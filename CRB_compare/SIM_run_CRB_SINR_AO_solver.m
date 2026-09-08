function [sol, hist] = SIM_run_CRB_SINR_AO_solver(params, chIn, alg, initState, gammaTarget)
%RUN_CRB_SINR_AO_SOLVER Alternating CRB-SINR covariance/phase optimizer.
%
% Outer iteration:
%   1) fixed theta: solve the transmit-covariance SDR and recover W;
%   2) accept the W block only when it is feasible and does not increase CRB;
%   3) use only Phi_L to build true covariance-domain SINR reserve;
%   4) update only Phi_1,...,Phi_{L-1} by feasible monotone CRB FP-SCA;
%   5) accept the complete phase block only when the final theta is feasible
%      and does not increase the original CRB.
%
% Every accepted phase trial is checked using the exact unregularized CRB
% and the exact SINR constraints. The AO layer then performs a second complete
% metric check before retaining the phase block.

ch = chIn;
ch.sigma_c2 = initState.sigma_c2;
ch.sigma_s2 = initState.sigma_s2;

theta = initState.theta0;
S = initState.S;
gammaTarget = gammaTarget(:);

maxOuter = alg.maxOuter;
blockAcceptRelTol = get_option(alg,'blockAcceptRelTol',1e-10);
maxSDRRetry = get_option(alg,'maxSDRRetry',2);
maxConsecutiveSDRFailures = get_option(alg,'maxConsecutiveSDRFailures',3);
progressEvery = get_option(alg,'progressEvery',100);
verbose = isfield(alg,'verbose') && logical(alg.verbose);
useOuterStop = logical(get_option(alg,'useOuterStop',true));
outerTol = get_option(alg,'tolOuter',1e-5);
outerMinIterations = max(1,round(get_option( ...
    alg,'outerMinIterations',15)));
outerStableIterations = max(1,round(get_option( ...
    alg,'outerStableIterations',6)));

hist.diagnosticEnabled = true;
hist.deltaTheta = NaN(maxOuter,1);
hist.deltaW = NaN(maxOuter,1);
hist.deltaP = NaN(maxOuter,1);
hist.branchSwitchFlag = false(maxOuter,1);
hist.WJumpRatio = NaN(maxOuter,1);
hist.PJumpRatio = NaN(maxOuter,1);
hist.ThetaJumpRatio = NaN(maxOuter,1);
hist.CRBBeforePhase = NaN(maxOuter,1);
hist.CRBAfterPhase = NaN(maxOuter,1);
hist.rateBeforePhase = NaN(maxOuter,1);
hist.rateAfterPhase = NaN(maxOuter,1);
hist.nmmseBeforePhase = NaN(maxOuter,1);
hist.nmmseAfterPhase = NaN(maxOuter,1);
hist.minSINRBeforePhase = NaN(maxOuter,1);
hist.minSINRAfterPhase = NaN(maxOuter,1);
hist.maxViolationBeforePhase = NaN(maxOuter,1);
hist.maxViolationAfterPhase = NaN(maxOuter,1);
hist.relMinEigABeforePhase = NaN(maxOuter,1);
hist.relMinEigBBeforePhase = NaN(maxOuter,1);
hist.relMinEigAAfterPhase = NaN(maxOuter,1);
hist.relMinEigBAfterPhase = NaN(maxOuter,1);
hist.traceAinvBeforePhase = NaN(maxOuter,1);
hist.traceBinvBeforePhase = NaN(maxOuter,1);
hist.traceAinvAfterPhase = NaN(maxOuter,1);
hist.traceBinvAfterPhase = NaN(maxOuter,1);
hist.conditionWarningBeforePhase = false(maxOuter,1);
hist.conditionWarningAfterPhase = false(maxOuter,1);
hist.phaseConditionRejects = zeros(maxOuter,1);
hist.phaseAbsoluteConditionRejects = zeros(maxOuter,1);
hist.phaseRelativeConditionRejects = zeros(maxOuter,1);
hist.phaseGuardTriggered = false(maxOuter,1);
hist.phaseMaxDropA = ones(maxOuter,1);
hist.phaseMaxDropB = ones(maxOuter,1);
hist.power = NaN(maxOuter,1);
hist.phaseInner = NaN(maxOuter,1);
hist.phasePenalty = NaN(maxOuter,1); % compatibility field; unused by FP-SCA
hist.phaseGradNorm = NaN(maxOuter,1);
hist.phaseAlpha = NaN(maxOuter,1);
hist.phaseSubproblemSolves = zeros(maxOuter,1);
hist.phaseModelResolveCount = zeros(maxOuter,1);
hist.phaseLineSearchTrials = zeros(maxOuter,1);
hist.phaseModelRejects = zeros(maxOuter,1);
hist.phaseFeasibilityRejects = zeros(maxOuter,1);
hist.phaseObjectiveRejects = zeros(maxOuter,1);
hist.phaseTrustRadiusEnd = NaN(maxOuter,1);
hist.phaseCurvatureFEnd = NaN(maxOuter,1);
hist.phaseCurvatureCMaxEnd = NaN(maxOuter,1);
hist.phasePredictedDecrease = NaN(maxOuter,1);
hist.phaseActualDecrease = NaN(maxOuter,1);
hist.outerRelativeDecrease = NaN(maxOuter,1);
hist.outerStableCount = zeros(maxOuter,1);
hist.WBlockAccepted = false(maxOuter,1);
hist.phaseBlockAccepted = false(maxOuter,1);
hist.WUpdateAttempted = false(maxOuter,1);
hist.WUpdateAccepted = false(maxOuter,1);
hist.PhaseUpdateAccepted = false(maxOuter,1);
hist.sdrFallbackUsed = false(maxOuter,1);
hist.sdrRetryCount = zeros(maxOuter,1);
hist.consecutiveSDRFailures = zeros(maxOuter,1);
hist.sdrStatus = cell(maxOuter,1);
hist.sdrCVXStatus = cell(maxOuter,1);
hist.recoveryValid = false(maxOuter,1);
hist.recoveryStatus = cell(maxOuter,1);
hist.recoveryCovarianceGapRel = NaN(maxOuter,1);
hist.recoveryPowerGapRel = NaN(maxOuter,1);
hist.recoveryDesiredGapRel = NaN(maxOuter,1);
hist.recoverySelfLeakRel = NaN(maxOuter,1);
hist.recoveryMinResidualRelEig = NaN(maxOuter,1);
hist.phaseStatus = cell(maxOuter,1);
hist.outerMarginAcceptedSteps = zeros(maxOuter,1);
hist.outerMarginStart = NaN(maxOuter,1);
hist.outerMarginEnd = NaN(maxOuter,1);
hist.outerMarginStatus = cell(maxOuter,1);
hist.outerMarginCRBDrift = NaN(maxOuter,1);

used = 0;
status = 'Not started';
failureStage = 'none';
failedOuterIteration = NaN;
firstSDRFeasible = false;
firstSDRStatus = 'not-run';
hadLaterSDRFailure = false;
consecutiveSDRFailures = 0;

bestFeasible = [];
bestCRB = Inf;
currentAccepted = [];
previousOuterCRB = Inf;
outerStableCount = 0;
outerStopTriggered = false;

% Before the first fixed-phase SDR at a new SINR target, use the inherited
% covariance only to retune the CRB-invariant outermost phase. This can turn
% a tight previous-target state into a useful communication carrier without
% any rank-one factorization or sensing-objective change.
hist.preOuterMarginApplied = false;
hist.preOuterMarginStart = NaN;
hist.preOuterMarginEnd = NaN;
hist.preOuterMarginStatus = 'not-run';
if isfield(initState,'Qcomm0') && ~isempty(initState.Qcomm0) && ...
        isfield(initState,'Rs0') && ~isempty(initState.Rs0) && ...
        logical(get_option(alg,'outerMarginEnable',true))
    [thetaPrepared,preOuterInfo] = SIM_update_outer_phase_margin( ...
        params,ch,theta,initState.Qcomm0,initState.Rs0,gammaTarget,alg);
    hist.preOuterMarginStart = get_option(preOuterInfo,'marginStart',NaN);
    hist.preOuterMarginEnd = get_option(preOuterInfo,'marginEnd',NaN);
    hist.preOuterMarginStatus = status_text(preOuterInfo,'status');
    if get_option(preOuterInfo,'acceptedSteps',0) > 0
        theta = thetaPrepared;
        initState.theta0 = thetaPrepared;
        hist.preOuterMarginApplied = true;
    end
end

% Cross-target continuation safeguard. When initState comes from the
% previous formal SINR point, evaluate its complete W/theta pair directly at
% the new target. If it is still feasible, use it as the initial incumbent.
% The first fixed-theta SDR solution may replace it only when the recovered
% complete solution does not increase the exact unregularized CRB.
[initialIncumbent,initialInfo] = build_initial_incumbent( ...
    params,ch,alg,initState,S,gammaTarget);
hist.initialIncumbentAvailable = initialInfo.available;
hist.initialIncumbentFeasibleAtTarget = initialInfo.feasible;
hist.initialIncumbentUsed = false;
hist.initialIncumbentCRB = initialInfo.CRB;
hist.initialIncumbentMinSINR = initialInfo.minSINR;
hist.initialIncumbentMaxViolation = initialInfo.maxViolation;
hist.initialIncumbentStatus = initialInfo.status;

if initialInfo.feasible
    currentAccepted = initialIncumbent;
    bestFeasible = initialIncumbent;
    bestCRB = initialIncumbent.metrics.CRB;
    previousOuterCRB = bestCRB;
    hist.initialIncumbentUsed = true;
end

for it = 1:maxOuter
    [P,~] = build_P(theta,ch.Omega);

    [sdr,sdrRetryCount] = solve_sdr_with_retry( ...
        params,ch,P,gammaTarget,alg,maxSDRRetry);
    hist.sdrRetryCount(it) = sdrRetryCount;
    hist.sdrStatus{it} = status_text(sdr,'status');
    hist.sdrCVXStatus{it} = status_text(sdr,'cvxStatus');

    recTrial = [];
    metricsTrial = [];
    sdrCandidateValid = false;

    if sdr.success
        recTrial = make_covariance_representation(sdr);
        hist.recoveryValid(it) = is_recovery_valid(recTrial);
        hist.recoveryStatus{it} = status_text(recTrial,'status');
        hist.recoveryCovarianceGapRel(it) = get_option( ...
            recTrial,'covarianceGapRel',NaN);
        hist.recoveryPowerGapRel(it) = get_option( ...
            recTrial,'powerGapRel',NaN);
        hist.recoveryDesiredGapRel(it) = get_option( ...
            recTrial,'maxDesiredSignalGapRel',NaN);
        hist.recoverySelfLeakRel(it) = get_option( ...
            recTrial,'maxSelfResidualLeakRel',NaN);
        hist.recoveryMinResidualRelEig(it) = get_option( ...
            recTrial,'minResidualRelEigBefore',NaN);

        if is_recovery_valid(recTrial)
            metricsTrial = SIM_evaluate_CRB_SINR_metrics( ...
                params,ch,P,sdr.Qcomm,sdr.Rs,S,recTrial.RxRecovered, ...
                gammaTarget,alg);
            sdrCandidateValid = is_feasible_metrics(metricsTrial);
        end
    end

    if it == 1
        firstSDRFeasible = sdr.success && sdrCandidateValid;
        if sdr.success && ~is_recovery_valid(recTrial)
            firstSDRStatus = ['Recovery invalid after ', ...
                status_text(sdr,'status'),': ',status_text(recTrial,'status')];
        elseif sdr.success && ~sdrCandidateValid
            firstSDRStatus = ['Recovered metrics invalid after ', ...
                status_text(sdr,'status')];
        else
            firstSDRStatus = status_text(sdr,'status');
        end
    end

    if sdr.success && sdrCandidateValid
        consecutiveSDRFailures = 0;
        wCandidate = make_candidate(theta,P,sdr,recTrial, ...
            metricsTrial,it,'before-phase');

        hist.WUpdateAttempted(it) = true;
        if isempty(currentAccepted) || is_not_worse( ...
                metricsTrial.CRB,currentAccepted.metrics.CRB, ...
                blockAcceptRelTol)
            if ~isempty(currentAccepted)
                oldRx = currentAccepted.rec.RxRecovered;
                hist.deltaW(it) = norm(wCandidate.rec.RxRecovered-oldRx,'fro');
                hist.WJumpRatio(it) = hist.deltaW(it)/max(norm(oldRx,'fro'),eps);
            else
                hist.deltaW(it) = NaN;
            end
            currentAccepted = wCandidate;
            hist.WBlockAccepted(it) = true;
            hist.WUpdateAccepted(it) = true;
        else
            % The old W/Rx pair is feasible for the same theta. Retain it
            % when the newly recovered SDR solution gives a larger CRB.
            hist.WBlockAccepted(it) = false;
        end
    else
        consecutiveSDRFailures = consecutiveSDRFailures + 1;
        hadLaterSDRFailure = hadLaterSDRFailure || it > 1;

        if isempty(currentAccepted)
            if it == 1
                failureStage = 'first-SDR';
            else
                failureStage = 'later-SDR';
            end
            failedOuterIteration = it;
            if sdr.success && ~is_recovery_valid(recTrial)
                status = sprintf(['Recovery was invalid at outer iteration ', ...
                    '%d after CVX status %s: %s'],it, ...
                    status_text(sdr,'status'),status_text(recTrial,'status'));
            elseif sdr.success
                status = sprintf(['Recovered metrics were invalid at outer ', ...
                    'iteration %d after CVX status %s.'], ...
                    it,status_text(sdr,'status'));
            else
                status = sprintf('SDR failed at outer iteration %d: %s', ...
                    it,status_text(sdr,'status'));
            end
            break;
        end

        hist.sdrFallbackUsed(it) = true;
        if consecutiveSDRFailures >= maxConsecutiveSDRFailures
            failureStage = 'later-SDR';
            failedOuterIteration = it;
            status = sprintf(['SDR failed for %d consecutive outer ', ...
                'iterations; retained the best earlier feasible pair. ', ...
                'Last status: %s'],consecutiveSDRFailures, ...
                status_text(sdr,'status'));
            break;
        end
        % Continue the phase block using the already accepted feasible W,
        % Rx, and theta rather than terminating on one numerical SDR failure.
    end

    hist.consecutiveSDRFailures(it) = consecutiveSDRFailures;

    metricsBefore = currentAccepted.metrics;
    Qcomm = currentAccepted.sdr.Qcomm;
    Rs = currentAccepted.sdr.Rs;
    Rx = currentAccepted.rec.RxRecovered;
    P = currentAccepted.P;

    used = it;
    hist.CRBBeforePhase(it) = metricsBefore.CRB;
    hist.rateBeforePhase(it) = metricsBefore.rate;
    hist.nmmseBeforePhase(it) = metricsBefore.nmmse;
    hist.minSINRBeforePhase(it) = min(metricsBefore.sinr);
    hist.maxViolationBeforePhase(it) = ...
        metricsBefore.maxScaledConstraintViolation;
    hist.relMinEigABeforePhase(it) = metricsBefore.crbInfo.relMinEigA;
    hist.relMinEigBBeforePhase(it) = metricsBefore.crbInfo.relMinEigB;
    hist.traceAinvBeforePhase(it) = metricsBefore.crbInfo.traceAinv;
    hist.traceBinvBeforePhase(it) = metricsBefore.crbInfo.traceBinv;
    hist.conditionWarningBeforePhase(it) = ...
        logical(metricsBefore.crbInfo.conditionWarning);
    hist.power(it) = real(trace(Rx));

    thetaOld = theta;
    POldForDiag = P;

    [thetaOuter,outerMarginInfo] = SIM_update_outer_phase_margin( ...
        params,ch,thetaOld,Qcomm,Rs,gammaTarget,alg);
    [Pouter,~] = build_P(thetaOuter,ch.Omega);
    metricsOuterTrial = SIM_evaluate_CRB_SINR_metrics( ...
        params,ch,Pouter,Qcomm,Rs,S,Rx,gammaTarget,alg);
    outerMarginAccepted = ...
        get_option(outerMarginInfo,'acceptedSteps',0) > 0 && ...
        is_feasible_metrics(metricsOuterTrial) && ...
        is_not_worse(metricsOuterTrial.CRB,metricsBefore.CRB, ...
        get_option(alg,'outerMarginCRBRelTol',1e-7));
    hist.outerMarginAcceptedSteps(it) = get_option( ...
        outerMarginInfo,'acceptedSteps',0);
    hist.outerMarginStart(it) = get_option(outerMarginInfo,'marginStart',NaN);
    hist.outerMarginEnd(it) = get_option(outerMarginInfo,'marginEnd',NaN);
    hist.outerMarginStatus{it} = status_text(outerMarginInfo,'status');
    hist.outerMarginCRBDrift(it) = get_option( ...
        outerMarginInfo,'crbRelativeDriftMax',NaN);

    if outerMarginAccepted
        thetaMiddle = thetaOuter;
        PMiddle = Pouter;
        metricsMiddle = metricsOuterTrial;
        currentAccepted = make_candidate(thetaMiddle,PMiddle, ...
            currentAccepted.sdr,currentAccepted.rec,metricsMiddle, ...
            it,'after-outer-margin');
    else
        thetaMiddle = thetaOld;
        PMiddle = P;
        metricsMiddle = metricsBefore;
    end

    [thetaInner,phaseInfo] = SIM_update_CRB_phase_FP_SCA( ...
        params,ch,thetaMiddle,Qcomm,Rs,gammaTarget,alg);
    [Pinner,~] = build_P(thetaInner,ch.Omega);
    metricsInnerTrial = SIM_evaluate_CRB_SINR_metrics( ...
        params,ch,Pinner,Qcomm,Rs,S,Rx,gammaTarget,alg);
    innerPhaseAccepted = get_option(phaseInfo,'acceptedSteps',0) > 0 && ...
        is_feasible_metrics(metricsInnerTrial) && ...
        is_not_worse(metricsInnerTrial.CRB,metricsMiddle.CRB, ...
        blockAcceptRelTol);

    if innerPhaseAccepted
        thetaNew = thetaInner;
        Pnew = Pinner;
        metricsAfter = metricsInnerTrial;
        currentAccepted = make_candidate(thetaNew,Pnew, ...
            currentAccepted.sdr,currentAccepted.rec,metricsAfter, ...
            it,'after-inner-phase');
    else
        thetaNew = thetaMiddle;
        Pnew = PMiddle;
        metricsAfter = metricsMiddle;
    end
    phaseAccepted = outerMarginAccepted || innerPhaseAccepted;

    if phaseAccepted
        hist.deltaTheta(it) = norm(exp(1j*thetaNew)-exp(1j*thetaOld),'fro');
        hist.ThetaJumpRatio(it) = hist.deltaTheta(it)/max(norm(exp(1j*thetaOld),'fro'),eps);
        hist.deltaW(it) = 0;
        hist.WJumpRatio(it) = 0;
        hist.deltaP(it) = norm(Pnew-POldForDiag,'fro');
        hist.PJumpRatio(it) = hist.deltaP(it)/max(norm(POldForDiag,'fro'),eps);
        hist.branchSwitchFlag(it) = hist.ThetaJumpRatio(it) > 0.5 || hist.WJumpRatio(it) > 0.5;
        theta = thetaNew;
        hist.phaseBlockAccepted(it) = true;
        hist.PhaseUpdateAccepted(it) = true;

        if is_better_return_candidate(currentAccepted,bestFeasible, ...
                blockAcceptRelTol)
            bestCRB = metricsAfter.CRB;
            bestFeasible = currentAccepted;
        end
    else
        theta = thetaOld;
        currentAccepted = make_candidate(thetaOld,P, ...
            currentAccepted.sdr,currentAccepted.rec,metricsBefore, ...
            it,'before-phase');
        metricsAfter = metricsBefore;
        hist.phaseBlockAccepted(it) = false;
    end

    if hist.WBlockAccepted(it) && metricsBefore.CRB < bestCRB
        bestCRB = metricsBefore.CRB;
        bestFeasible = make_candidate(thetaOld,P, ...
            currentAccepted.sdr,currentAccepted.rec,metricsBefore, ...
            it,'before-phase');
    end

    % A successful phase move changes theta and gives the next SDR a new
    % numerical problem, so prior consecutive SDR failures no longer count.
    if phaseAccepted
        consecutiveSDRFailures = 0;
    end

    hist.CRBAfterPhase(it) = metricsAfter.CRB;
    hist.rateAfterPhase(it) = metricsAfter.rate;
    hist.nmmseAfterPhase(it) = metricsAfter.nmmse;
    hist.minSINRAfterPhase(it) = min(metricsAfter.sinr);
    hist.maxViolationAfterPhase(it) = ...
        metricsAfter.maxScaledConstraintViolation;
    hist.relMinEigAAfterPhase(it) = metricsAfter.crbInfo.relMinEigA;
    hist.relMinEigBAfterPhase(it) = metricsAfter.crbInfo.relMinEigB;
    hist.traceAinvAfterPhase(it) = metricsAfter.crbInfo.traceAinv;
    hist.traceBinvAfterPhase(it) = metricsAfter.crbInfo.traceBinv;
    hist.conditionWarningAfterPhase(it) = ...
        logical(metricsAfter.crbInfo.conditionWarning);
    hist.phaseConditionRejects(it) = get_option( ...
        phaseInfo,'conditionRejects',0);
    hist.phaseAbsoluteConditionRejects(it) = get_option( ...
        phaseInfo,'absoluteConditionRejects',0);
    hist.phaseRelativeConditionRejects(it) = get_option( ...
        phaseInfo,'relativeConditionRejects',0);
    hist.phaseGuardTriggered(it) = logical(get_option( ...
        phaseInfo,'phaseGuardTriggered',false));
    hist.phaseMaxDropA(it) = get_option(phaseInfo,'maxPhaseDropA',1);
    hist.phaseMaxDropB(it) = get_option(phaseInfo,'maxPhaseDropB',1);
    hist.phaseInner(it) = get_option(phaseInfo,'innerUsed',0);
    hist.phasePenalty(it) = NaN;
    hist.phaseGradNorm(it) = get_option(phaseInfo,'lastGradNorm',NaN);
    hist.phaseAlpha(it) = get_option(phaseInfo,'lastAlpha',NaN);
    hist.phaseSubproblemSolves(it) = get_option( ...
        phaseInfo,'subproblemSolves',0);
    hist.phaseModelResolveCount(it) = get_option( ...
        phaseInfo,'modelResolveCount',0);
    hist.phaseLineSearchTrials(it) = get_option( ...
        phaseInfo,'lineSearchTrials',0);
    hist.phaseModelRejects(it) = get_option(phaseInfo,'modelRejects',0);
    hist.phaseFeasibilityRejects(it) = get_option( ...
        phaseInfo,'feasibilityRejects',0);
    hist.phaseObjectiveRejects(it) = get_option( ...
        phaseInfo,'objectiveRejects',0);
    hist.phaseTrustRadiusEnd(it) = get_option( ...
        phaseInfo,'trustRadiusEnd',NaN);
    hist.phaseCurvatureFEnd(it) = get_option( ...
        phaseInfo,'curvatureFEnd',NaN);
    hist.phaseCurvatureCMaxEnd(it) = get_option( ...
        phaseInfo,'curvatureCMaxEnd',NaN);
    hist.phasePredictedDecrease(it) = get_option( ...
        phaseInfo,'lastPredictedDecrease',NaN);
    hist.phaseActualDecrease(it) = get_option( ...
        phaseInfo,'lastActualDecrease',NaN);
    hist.phaseStatus{it} = status_text(phaseInfo,'status');

    if isfinite(previousOuterCRB)
        hist.outerRelativeDecrease(it) = max(0, ...
            (previousOuterCRB-metricsAfter.CRB) / ...
            max(abs(previousOuterCRB),realmin));
        if hist.outerRelativeDecrease(it) <= outerTol
            outerStableCount = outerStableCount + 1;
        else
            outerStableCount = 0;
        end
    else
        outerStableCount = 0;
    end
    hist.outerStableCount(it) = outerStableCount;
    previousOuterCRB = metricsAfter.CRB;

    if verbose || mod(it,progressEvery)==0 || it==1 || it==maxOuter
        fprintf(['    outer %04d/%04d: CRB %.3e -> %.3e, ', ...
            'minSINR %.3e -> %.3e, relEig(A/B)=%.1e/%.1e, ', ...
            'Wacc=%d, outerMargin=%d, innerPhi=%d, retry=%d, guardReject=%d\n'], ...
            it,maxOuter,metricsBefore.CRB,metricsAfter.CRB, ...
            min(metricsBefore.sinr),min(metricsAfter.sinr), ...
            metricsAfter.crbInfo.relMinEigA, ...
            metricsAfter.crbInfo.relMinEigB, ...
            hist.WBlockAccepted(it),outerMarginAccepted,innerPhaseAccepted, ...
            hist.sdrRetryCount(it), ...
            hist.phaseConditionRejects(it));
    end

    if useOuterStop && it >= outerMinIterations && ...
            outerStableCount >= outerStableIterations
        status = sprintf([ ...
            'Outer monotone convergence reached after %d iterations ', ...
            '(%d consecutive relative changes <= %.3e).'], ...
            it,outerStableCount,outerTol);
        outerStopTriggered = true;
        break;
    end
end

% Final fixed-theta consistency solve. It is another candidate only; a
% failure or a worse recovered solution does not discard an accepted pair.
[Pfinal,~] = build_P(theta,ch.Omega);
[sdrFinal,finalRetryCount] = solve_sdr_with_retry( ...
    params,ch,Pfinal,gammaTarget,alg,maxSDRRetry);

if sdrFinal.success
    recFinal = make_covariance_representation(sdrFinal);
    if is_recovery_valid(recFinal)
        metricsFinal = SIM_evaluate_CRB_SINR_metrics(params,ch,Pfinal, ...
            sdrFinal.Qcomm,sdrFinal.Rs,S,recFinal.RxRecovered,gammaTarget,alg);

        if is_feasible_metrics(metricsFinal) && metricsFinal.CRB < bestCRB
            bestCRB = metricsFinal.CRB;
            bestFeasible = make_candidate(theta,Pfinal,sdrFinal,recFinal, ...
                metricsFinal,used+1,'final-consistency');
        end
    else
        metricsFinal = [];
    end
else
    metricsFinal = [];
end

if isempty(bestFeasible) && ~isempty(currentAccepted)
    bestFeasible = currentAccepted;
    bestCRB = currentAccepted.metrics.CRB; %#ok<NASGU>
end

if isempty(bestFeasible)
    if strcmp(failureStage,'none')
        failureStage = 'final-SDR';
    end
    if isempty(status) || strcmp(status,'Not started')
        status = sprintf('No feasible CRB-SINR pair was retained. Final SDR: %s', ...
            status_text(sdrFinal,'status'));
    end
else
    if strcmp(status,'Not started')
        if hadLaterSDRFailure
            status = ['Solved with later SDR fallbacks; returned the best ', ...
                'feasible candidate.'];
        elseif sdrFinal.success
            status = 'Solved; returned the best feasible candidate.';
        else
            status = ['Final consistency SDR failed; returned the best ', ...
                'earlier feasible candidate.'];
        end
    end
end

hist = trim_history(hist,used);
hist.outerIterationsCompleted = used;
hist.bestIteration = NaN;
hist.bestStage = 'none';
hist.finalSDRStatus = status_text(sdrFinal,'status');
hist.finalSDRRetryCount = finalRetryCount;
hist.failureStage = failureStage;
hist.failedOuterIteration = failedOuterIteration;
hist.firstSDRFeasible = firstSDRFeasible;
hist.firstSDRStatus = firstSDRStatus;
hist.hadLaterSDRFailure = hadLaterSDRFailure;
hist.outerStopTriggered = outerStopTriggered;

if isempty(bestFeasible)
    sol = make_failure_solution(params,theta,S,ch,gammaTarget, ...
        status,failureStage,failedOuterIteration,firstSDRFeasible, ...
        firstSDRStatus);
    return;
end

sol.success = true;
sol.status = status;
% Always return the best complete feasible pair found over all AO stages,
% not the final iterate. The next SINR target therefore inherits the phase
% matrix that belongs to the lowest exact physical CRB retained at this
% target, together with its matched beamformer/covariance solution.
sol.W = [];
sol.Wc = [];
sol.Wr = [];
sol.theta = bestFeasible.theta;
sol.P = bestFeasible.P;
sol.Rx = bestFeasible.rec.RxRecovered;
sol.bestW = [];
sol.bestTheta = sol.theta;
sol.bestP = sol.P;
sol.bestRx = sol.Rx;
sol.QcommSDR = bestFeasible.sdr.Qcomm;
sol.RsSDR = bestFeasible.sdr.Rs;
sol.RsRecovered = [];
sol.recoveryValid = is_recovery_valid(bestFeasible.rec);
sol.recoveryStatus = status_text(bestFeasible.rec,'status');
sol.recoveryCovarianceGapRel = get_option( ...
    bestFeasible.rec,'covarianceGapRel',NaN);
sol.recoveryPowerGapRel = get_option(bestFeasible.rec,'powerGapRel',NaN);
sol.recoveryDesiredGapRel = get_option( ...
    bestFeasible.rec,'maxDesiredSignalGapRel',NaN);
sol.recoverySelfLeakRel = get_option( ...
    bestFeasible.rec,'maxSelfResidualLeakRel',NaN);
sol.recoveryMinResidualRelEig = get_option( ...
    bestFeasible.rec,'minResidualRelEigBefore',NaN);
sol.S = S;
sol.sigma_c2 = ch.sigma_c2;
sol.sigma_s2 = ch.sigma_s2;
sol.metrics = bestFeasible.metrics;
sol.gammaTarget = gammaTarget;
sol.lambda = [];
sol.penalty = NaN;
sol.params = params;
sol.bestIteration = bestFeasible.iteration;
sol.bestStage = bestFeasible.stage;
sol.failureStage = failureStage;
sol.failedOuterIteration = failedOuterIteration;
sol.firstSDRFeasible = firstSDRFeasible;
sol.firstSDRStatus = firstSDRStatus;
sol.initialIncumbentAvailable = initialInfo.available;
sol.initialIncumbentFeasibleAtTarget = initialInfo.feasible;
sol.initialIncumbentUsed = hist.initialIncumbentUsed;
sol.initialIncumbentCRB = initialInfo.CRB;
sol.initialIncumbentMinSINR = initialInfo.minSINR;
sol.initialIncumbentMaxViolation = initialInfo.maxViolation;
sol.phaseGuardTriggered = any(hist.phaseGuardTriggered);
sol.phaseGuardRejectCount = sum(hist.phaseConditionRejects);
sol.phaseAbsoluteGuardRejectCount = ...
    sum(hist.phaseAbsoluteConditionRejects);
sol.phaseRelativeGuardRejectCount = ...
    sum(hist.phaseRelativeConditionRejects);
if isempty(hist.phaseMaxDropA)
    sol.phaseGuardMaxDropA = 1;
    sol.phaseGuardMaxDropB = 1;
else
    sol.phaseGuardMaxDropA = max(hist.phaseMaxDropA);
    sol.phaseGuardMaxDropB = max(hist.phaseMaxDropB);
end
sol.phaseGuardMinRelEigA = finite_min([ ...
    hist.relMinEigABeforePhase;hist.relMinEigAAfterPhase],NaN);
sol.phaseGuardMinRelEigB = finite_min([ ...
    hist.relMinEigBBeforePhase;hist.relMinEigBAfterPhase],NaN);

hist.bestIteration = bestFeasible.iteration;
hist.bestStage = bestFeasible.stage;
hist.finalMetrics = bestFeasible.metrics;
end

function [candidate,info] = build_initial_incumbent( ...
    params,ch,alg,initState,S,gammaTarget)
%BUILD_INITIAL_INCUMBENT Re-evaluate the complete warm-start W/theta pair.
%
% This function does not optimize anything. It only reconstructs the
% covariance representation associated with initState.W0 and evaluates the
% exact CRB/SINR metrics at the new target. Consequently, a feasible warm
% start provides a rigorous cross-target CRB incumbent.

candidate = [];
info.available = false;
info.feasible = false;
info.CRB = NaN;
info.minSINR = NaN;
info.maxViolation = NaN;
info.status = 'warm-start-complete-pair-unavailable';

if ~isstruct(initState) || ~isfield(initState,'Qcomm0') || ...
        isempty(initState.Qcomm0) || ~isfield(initState,'Rs0') || ...
        isempty(initState.Rs0) || ~isfield(initState,'theta0') || ...
        isempty(initState.theta0)
    return;
end
Qcomm = initState.Qcomm0;
Rs = initState.Rs0;
if size(Qcomm,1) ~= params.Nt || size(Qcomm,3) ~= params.K || ...
        ~isequal(size(Rs),[params.Nt,params.Nt])
    info.status = 'warm-start-covariances-have-incompatible-dimensions';
    return;
end

try
    [P,~] = build_P(initState.theta0,ch.Omega);
catch ME
    info.status = ['warm-start-phase-build-failed: ',ME.message];
    return;
end

sdr.success = true;
sdr.status = 'covariance-warm-start-incumbent';
sdr.cvxStatus = 'not-applicable';
sdr.Qcomm = Qcomm;
sdr.Rs = Rs;
rec = make_covariance_representation(sdr);
Rx = rec.RxRecovered;

try
    metrics = SIM_evaluate_CRB_SINR_metrics( ...
        params,ch,P,Qcomm,Rs,S,Rx,gammaTarget,alg);
catch ME
    info.available = true;
    info.status = ['warm-start-metric-evaluation-failed: ',ME.message];
    return;
end

info.available = true;
info.CRB = metrics.CRB;
if isfield(metrics,'sinr') && ~isempty(metrics.sinr)
    info.minSINR = min(metrics.sinr);
end
if isfield(metrics,'maxScaledConstraintViolation')
    info.maxViolation = metrics.maxScaledConstraintViolation;
end
info.feasible = is_feasible_metrics(metrics);
if info.feasible
    info.status = 'complete-warm-start-feasible-at-new-target';
    candidate = make_candidate(initState.theta0,P,sdr,rec,metrics, ...
        0,'initial-incumbent');
else
    info.status = 'complete-warm-start-infeasible-at-new-target';
end
end

function [sdr,retryCount] = solve_sdr_with_retry( ...
    params,ch,P,gammaTarget,alg,maxRetry)
retryCount = 0;
% Fixed-schema seed; deliberately avoid empty-structure assignment.
sdr = struct('success',false,'status','not-run','cvxStatus','not-run');
for ir = 0:maxRetry
    if ir > 0
        retryCount = ir;
    end
    try
        cvx_clear;
    catch
    end
    if logical(get_option(alg,'useRobustSDR',false))
        sdr = SIM_solve_CRB_SINR_SDR_robust( ...
            params,ch,P,gammaTarget,alg);
    else
        sdr = SIM_solve_CRB_SINR_SDR(params,ch,P,gammaTarget,alg);
    end
    if isfield(sdr,'success') && isscalar(sdr.success) && logical(sdr.success)
        return;
    end
end
end

function candidate = make_candidate(theta,P,sdr,rec,metrics,iteration,stage)
candidate.theta = theta;
candidate.P = P;
candidate.sdr = sdr;
candidate.rec = rec;
candidate.metrics = metrics;
candidate.iteration = iteration;
candidate.stage = stage;
end

function tf = is_not_worse(newValue,oldValue,relTol)
scale = max(abs(oldValue),realmin);
tf = isfinite(newValue) && isfinite(oldValue) && ...
    newValue <= oldValue + relTol*scale;
end

function tf = is_better_return_candidate(candidate,incumbent,relTol)
% Prefer lower CRB. When CRBs are numerically indistinguishable, retain the
% phase/covariance pair with the larger true minimum SINR so that the saved
% continuation state carries useful communication reserve.
tf = false;
if isempty(candidate) || ~isfield(candidate,'metrics')
    return;
end
if isempty(incumbent) || ~isfield(incumbent,'metrics')
    tf = true;
    return;
end
newCRB = candidate.metrics.CRB;
oldCRB = incumbent.metrics.CRB;
scale = max(abs(oldCRB),realmin);
tol = relTol*scale;
if newCRB < oldCRB-tol
    tf = true;
elseif abs(newCRB-oldCRB) <= tol
    tf = min(candidate.metrics.sinr) > min(incumbent.metrics.sinr);
end
end

function tf = is_recovery_valid(rec)
tf = isstruct(rec) && isfield(rec,'valid') && isscalar(rec.valid) && ...
    isfield(rec,'RxRecovered') && ~isempty(rec.RxRecovered) && ...
    logical(rec.valid) && all(isfinite(rec.RxRecovered(:)));
end

function rec = make_covariance_representation(sdr)
% Fixed-schema adapter used by the AO bookkeeping. It does not recover W.
rec.valid = isstruct(sdr) && isfield(sdr,'success') && logical(sdr.success);
rec.status = 'unrecovered-SDR-covariances';
rec.W = []; rec.Wc = []; rec.Wr = []; rec.RsRecovered = [];
rec.RxRecovered = sdr.Rs;
for kk=1:size(sdr.Qcomm,3), rec.RxRecovered=rec.RxRecovered+sdr.Qcomm(:,:,kk); end
rec.RxRecovered=(rec.RxRecovered+rec.RxRecovered')/2;
rec.RxSDR=rec.RxRecovered; rec.residuals=[];
rec.power=real(trace(rec.RxRecovered)); rec.powerSDR=rec.power;
rec.covarianceGapRel=0; rec.powerGapRel=0;
rec.maxDesiredSignalGapRel=0; rec.maxSelfResidualLeakRel=0;
rec.minResidualRelEigBefore=NaN;
end

function tf = is_feasible_metrics(metrics)
tf = isstruct(metrics) ...
    && isfield(metrics,'CRB') && isscalar(metrics.CRB) ...
    && isfinite(metrics.CRB) && metrics.CRB > 0 ...
    && isfield(metrics,'crbValid') && isscalar(metrics.crbValid) ...
    && logical(metrics.crbValid) ...
    && isfield(metrics,'sinrFeasible') && isscalar(metrics.sinrFeasible) ...
    && logical(metrics.sinrFeasible);
end

function sol = make_failure_solution(params,theta,S,ch,gammaTarget, ...
    status,failureStage,failedOuterIteration,firstSDRFeasible, ...
    firstSDRStatus)
sol.success = false;
sol.status = status;
sol.W = [];
sol.Wc = [];
sol.Wr = [];
sol.theta = theta;
sol.P = [];
sol.Rx = [];
sol.QcommSDR = [];
sol.RsSDR = [];
sol.RsRecovered = [];
sol.recoveryValid = false;
sol.recoveryStatus = 'no-valid-recovery';
sol.recoveryCovarianceGapRel = NaN;
sol.recoveryPowerGapRel = NaN;
sol.recoveryDesiredGapRel = NaN;
sol.recoverySelfLeakRel = NaN;
sol.recoveryMinResidualRelEig = NaN;
sol.S = S;
sol.sigma_c2 = ch.sigma_c2;
sol.sigma_s2 = ch.sigma_s2;
sol.metrics = struct('CRB',Inf,'crbValid',false,'sinrFeasible',false, ...
    'sinr',zeros(params.K,1),'crbInfo',struct('valid',false, ...
    'relMinEigA',NaN,'relMinEigB',NaN,'rcondA',NaN,'rcondB',NaN));
sol.gammaTarget = gammaTarget;
sol.lambda = [];
sol.penalty = NaN;
sol.params = params;
sol.bestIteration = NaN;
sol.bestStage = 'none';
sol.failureStage = failureStage;
sol.failedOuterIteration = failedOuterIteration;
sol.firstSDRFeasible = firstSDRFeasible;
sol.firstSDRStatus = firstSDRStatus;
sol.initialIncumbentAvailable = false;
sol.initialIncumbentFeasibleAtTarget = false;
sol.initialIncumbentUsed = false;
sol.initialIncumbentCRB = NaN;
sol.initialIncumbentMinSINR = NaN;
sol.initialIncumbentMaxViolation = NaN;
sol.phaseGuardTriggered = false;
sol.phaseGuardRejectCount = 0;
sol.phaseAbsoluteGuardRejectCount = 0;
sol.phaseRelativeGuardRejectCount = 0;
sol.phaseGuardMaxDropA = 1;
sol.phaseGuardMaxDropB = 1;
sol.phaseGuardMinRelEigA = NaN;
sol.phaseGuardMinRelEigB = NaN;
end

function value = finite_min(values,defaultValue)
values = values(isfinite(values));
if isempty(values)
    value = defaultValue;
else
    value = min(values);
end
end

function hist = trim_history(hist,used)
fields = fieldnames(hist);
for ii = 1:numel(fields)
    name = fields{ii};
    value = hist.(name);
    if iscell(value) && isvector(value) && numel(value) >= used
        hist.(name) = value(1:used);
    elseif (isnumeric(value) || islogical(value)) && isvector(value) ...
            && numel(value) >= used
        hist.(name) = value(1:used);
    end
end
hist.iteration = (1:used).';
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end

function textValue = status_text(s,fieldName)
textValue = 'not-reported';
if isstruct(s) && isfield(s,fieldName) && ~isempty(s.(fieldName))
    textValue = char(string(s.(fieldName)));
end
end
