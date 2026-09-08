function [theta,info] = SIM_update_CRB_phase_FP_SCA( ...
    params,ch,theta,Qcomm,Rs,gammaTarget,alg)
%SIM_UPDATE_CRB_PHASE_FP_SCA Feasible and monotone SIM phase update.
%
% For fixed covariance, this routine updates only the inner phase-angle
% blocks theta(:,1:L-1). The outermost layer Phi_L is deliberately excluded
% because it is CRB-invariant for the extended-target model and is handled
% separately by SIM_update_outer_phase_margin. Each
% accepted trial is checked using the exact unregularized CRB and the exact
% SINR inequalities. Therefore an accepted step satisfies both:
%
%   1) the CRB does not increase, within the configured numerical tolerance;
%   2) every user SINR constraint remains feasible, within the configured
%      numerical tolerance.
%
% A small convex proximal-SCA subproblem generates a candidate direction.
% When internal regularization is enabled, the local model and gradient use
% the same stabilized CRB surrogate, while acceptance still requires the
% exact unregularized CRB to be non-increasing. To control complexity, the
% direction is solved only a few times. Cheap exact line-search trials are
% attempted before the model is rebuilt.

Rx = Rs;
for kk = 1:params.K, Rx = Rx + Qcomm(:,:,kk); end
Rx = (Rx+Rx')/2;
gammaTarget = gammaTarget(:);
activeLayers = 1:max(0,ch.L-1);

maxInner = max(0,round(get_option(alg,'maxPhaseInner',1)));
maxModelResolve = max(1,round(get_option(alg,'phaseMaxModelResolve',2)));
maxLineSearch = max(1,round(get_option(alg,'phaseMaxLineSearch',8)));
lineSearchBeta = get_option(alg,'phaseLineSearchBeta',0.5);
minStepScale = get_option(alg,'phaseMinStepScale',1/128);
stepTol = get_option(alg,'phaseStepTol',1e-6);
gradTol = get_option(alg,'phaseGradTol',1e-7);
modelTol = get_option(alg,'phaseModelTol',1e-8);
objectiveTol = get_option(alg,'phaseObjectiveTol',1e-10);
sufficientDecrease = get_option(alg,'phaseSufficientDecrease',1e-4);
constraintTol = get_option(alg,'phaseExactConstraintTol', ...
    get_option(alg,'phaseConstraintTol',1e-5));
subproblemConstraintTol = get_option(alg, ...
    'phaseSubproblemConstraintTol',constraintTol);
requireModelUpperBound = logical(get_option( ...
    alg,'phaseRequireModelUpperBound',true));

trustRadius = get_option(alg,'phaseTrustRadius0',0.08);
trustRadiusMin = get_option(alg,'phaseTrustRadiusMin',1e-5);
trustRadiusMax = get_option(alg,'phaseTrustRadiusMax',0.20);
trustShrink = get_option(alg,'phaseTrustShrink',0.5);
trustExpand = get_option(alg,'phaseTrustExpand',1.10);

rhoF = get_option(alg,'phaseCurvatureF0',1);
rhoC = get_option(alg,'phaseCurvatureC0',1)*ones(params.K,1);
curvatureGrowth = get_option(alg,'phaseCurvatureGrowth',2);
curvatureDecrease = get_option(alg,'phaseCurvatureDecrease',1.1);
curvatureMin = get_option(alg,'phaseCurvatureMin',1e-4);
curvatureMax = get_option(alg,'phaseCurvatureMax',1e8);

phaseGuardTol = get_option(alg,'phaseGuardMinRelEig',0);
phaseMaxRelativeEigDrop = get_option( ...
    alg,'phaseMaxRelativeEigDrop',Inf);

[Pstart,~] = build_P(theta,ch.Omega);
[Jstart,~,crbInfoStart,searchLogStart] = SIM_crb_phase_value_gradient( ...
    params,ch,Pstart,theta,Rx,alg);
[~,~,baseConstraintInfo] = SIM_sinr_phase_constraints( ...
    params,ch,Pstart,theta,Qcomm,Rs,gammaTarget,[],false);
marginStart = min(baseConstraintInfo.sinr./max(gammaTarget,realmin)-1);
reserveFraction = get_option(alg,'innerReserveFraction',0.10);
reserveFloor = get_option(alg,'innerReserveFloor',1e-5);
retainedReserve = max(reserveFloor,reserveFraction*max(marginStart,0));
phaseGammaTarget = gammaTarget*(1+retainedReserve);
[cStart,~,constraintInfoStart] = SIM_sinr_phase_constraints( ...
    params,ch,Pstart,theta,Qcomm,Rs,phaseGammaTarget,[],false);
constraintScale = constraintInfoStart.scale;

info.success = crbInfoStart.valid;
info.status = 'not-started';
info.innerUsed = 0;
info.acceptedSteps = 0;
info.subproblemSolves = 0;
info.modelResolveCount = 0;
info.lineSearchTrials = 0;
info.modelRejects = 0;
info.feasibilityRejects = 0;
info.objectiveRejects = 0;
info.conditionRejects = 0;
info.absoluteConditionRejects = 0;
info.relativeConditionRejects = 0;
info.phaseGuardTriggered = false;
info.maxPhaseDropA = 1;
info.maxPhaseDropB = 1;
info.CRBStart = Jstart;
if ~isfinite(searchLogStart)
    info.success = false;
end
info.CRBEnd = Jstart;
info.maxViolationStart = max([0;cStart]);
info.maxViolationEnd = info.maxViolationStart;
info.lastAlpha = NaN;
info.lastGradNorm = NaN;
info.lastPredictedDecrease = 0;
info.lastActualDecrease = 0;
info.trustRadiusStart = trustRadius;
info.trustRadiusEnd = trustRadius;
info.curvatureFStart = rhoF;
info.curvatureFEnd = rhoF;
info.curvatureCMaxStart = max(rhoC);
info.curvatureCMaxEnd = max(rhoC);
info.relMinEigAStart = crbInfoStart.relMinEigA;
info.relMinEigBStart = crbInfoStart.relMinEigB;
info.relMinEigAEnd = crbInfoStart.relMinEigA;
info.relMinEigBEnd = crbInfoStart.relMinEigB;
info.activeLayers = activeLayers;
info.entryRelativeSINRMargin = marginStart;
info.retainedRelativeSINRReserve = retainedReserve;

if ~crbInfoStart.valid
    info.status = ['invalid initial CRB: ',crbInfoStart.invalidReason];
    return;
end
if info.maxViolationStart > constraintTol
    info.status = sprintf([ ...
        'initial phase is outside the phase feasibility tolerance: ', ...
        '%.3e > %.3e; phase retained.'], ...
        info.maxViolationStart,constraintTol);
    return;
end
if maxInner == 0
    info.status = 'phase update disabled by maxPhaseInner=0.';
    return;
end
if isempty(activeLayers)
    info.status = 'no inner sensing layers; CRB phase update skipped.';
    return;
end

for inner = 1:maxInner
    [Pcur,~] = build_P(theta,ch.Omega);
    [Jcur,~,crbInfoCur,searchLogCur] = SIM_crb_phase_value_gradient( ...
        params,ch,Pcur,theta,Rx,alg);
    if ~crbInfoCur.valid
        info.success = false;
        info.status = ['CRB became invalid: ',crbInfoCur.invalidReason];
        break;
    end

    gradFFull = crbInfoCur.gradLogTheta;
    gradF = gradFFull(:,activeLayers);
    fSearchCur = searchLogCur;
    fRawCur = crbInfoCur.logCRB;
    [cCur,gradCFull] = SIM_sinr_phase_constraints( ...
        params,ch,Pcur,theta,Qcomm,Rs,phaseGammaTarget,constraintScale,true);
    gradC = gradCFull(:,activeLayers,:);

    gradNorm = norm(gradF(:));
    info.lastGradNorm = gradNorm;
    if gradNorm <= gradTol
        info.status = 'phase log-CRB gradient tolerance reached.';
        break;
    end

    acceptedThisInner = false;
    noUsefulDirection = false;

    for modelResolve = 1:maxModelResolve
        subproblem = SIM_solve_phase_FP_SCA_subproblem( ...
            gradF,cCur,gradC,rhoF,rhoC,trustRadius, ...
            subproblemConstraintTol,alg);
        info.subproblemSolves = info.subproblemSolves + 1;
        if modelResolve > 1
            info.modelResolveCount = info.modelResolveCount + 1;
        end

        if ~subproblem.success
            rhoF = min(curvatureMax,curvatureGrowth*rhoF);
            rhoC = min(curvatureMax,curvatureGrowth*rhoC);
            trustRadius = max(trustRadiusMin,trustShrink*trustRadius);
            continue;
        end

        dTheta = subproblem.dTheta;
        if subproblem.stepInfNorm <= stepTol || ...
                subproblem.predictedDecrease <= objectiveTol
            noUsefulDirection = true;
            break;
        end

        alpha = 1;
        for lineTrial = 1:maxLineSearch
            info.lineSearchTrials = info.lineSearchTrials + 1;
            if alpha < minStepScale
                break;
            end

            step = alpha*dTheta;
            thetaTry = theta;
            thetaTry(:,activeLayers) = mod( ...
                theta(:,activeLayers)+step,2*pi);
            [Ptry,~] = build_P(thetaTry,ch.Omega);
            [Jtry,crbTryInfo] = SIM_crb_value(ch,Ptry,Rx,alg);
            searchLogTry = crbTryInfo.logCRB;
            if logical(get_option(alg,'useInternalCRBRegularization',false))
                [~,crbTryRegInfo] = SIM_crb_value_regularized( ...
                    ch,Ptry,Rx,alg);
                if crbTryRegInfo.valid
                    searchLogTry = crbTryRegInfo.logCRB;
                end
            end

            absoluteConditionOK = crbTryInfo.valid;
            if absoluteConditionOK && phaseGuardTol > 0
                dynamicGuardA = min(phaseGuardTol, ...
                    max(0,(1-1e-8)*crbInfoCur.relMinEigA));
                dynamicGuardB = min(phaseGuardTol, ...
                    max(0,(1-1e-8)*crbInfoCur.relMinEigB));
                absoluteConditionOK = ...
                    crbTryInfo.relMinEigA >= dynamicGuardA && ...
                    crbTryInfo.relMinEigB >= dynamicGuardB;
            end

            if ~absoluteConditionOK
                info.conditionRejects = info.conditionRejects + 1;
                info.absoluteConditionRejects = ...
                    info.absoluteConditionRejects + 1;
                info.phaseGuardTriggered = true;
                alpha = lineSearchBeta*alpha;
                continue;
            end

            phaseDropA = crbInfoCur.relMinEigA / ...
                max(crbTryInfo.relMinEigA,realmin);
            phaseDropB = crbInfoCur.relMinEigB / ...
                max(crbTryInfo.relMinEigB,realmin);
            info.maxPhaseDropA = max(info.maxPhaseDropA,phaseDropA);
            info.maxPhaseDropB = max(info.maxPhaseDropB,phaseDropB);
            relativeConditionOK = ...
                phaseDropA <= phaseMaxRelativeEigDrop && ...
                phaseDropB <= phaseMaxRelativeEigDrop;
            if ~relativeConditionOK
                info.conditionRejects = info.conditionRejects + 1;
                info.relativeConditionRejects = ...
                    info.relativeConditionRejects + 1;
                info.phaseGuardTriggered = true;
                alpha = lineSearchBeta*alpha;
                continue;
            end

            [cTry,~,~] = SIM_sinr_phase_constraints( ...
                params,ch,Ptry,thetaTry,Qcomm,Rs,phaseGammaTarget, ...
                constraintScale,false);
            exactFeasible = all(cTry <= constraintTol);
            if ~exactFeasible
                info.feasibilityRejects = info.feasibilityRejects + 1;
                alpha = lineSearchBeta*alpha;
                continue;
            end

            stepNormSquared = sum(step(:).^2);
            objectiveModel = fSearchCur + sum(gradF(:).*step(:)) + ...
                0.5*rhoF*stepNormSquared;
            constraintModel = zeros(params.K,1);
            badConstraintModel = false(params.K,1);
            for k = 1:params.K
                gradK = gradC(:,:,k);
                constraintModel(k) = cCur(k) + ...
                    sum(gradK(:).*step(:)) + ...
                    0.5*rhoC(k)*stepNormSquared;
                badConstraintModel(k) = ...
                    cTry(k) > constraintModel(k)+modelTol;
            end

            objectiveModelValid = ...
                searchLogTry <= objectiveModel+modelTol;
            constraintModelValid = ~any(badConstraintModel);
            modelValid = objectiveModelValid && constraintModelValid;

            modelChange = objectiveModel-fSearchCur;
            sufficientDescentOK = searchLogTry <= ...
                fSearchCur+sufficientDecrease*modelChange+objectiveTol;
            % The reported physical CRB remains the exact unregularized
            % objective and is never allowed to increase beyond tolerance.
            monotoneOK = crbTryInfo.logCRB <= fRawCur+objectiveTol;

            if requireModelUpperBound && ~modelValid
                info.modelRejects = info.modelRejects + 1;
                alpha = lineSearchBeta*alpha;
                continue;
            end
            if ~(sufficientDescentOK && monotoneOK)
                info.objectiveRejects = info.objectiveRejects + 1;
                alpha = lineSearchBeta*alpha;
                continue;
            end

            theta = thetaTry;
            info.acceptedSteps = info.acceptedSteps + 1;
            info.innerUsed = inner;
            info.lastAlpha = alpha;
            info.lastPredictedDecrease = max(0,-modelChange);
            info.lastActualDecrease = max(0,fRawCur-crbTryInfo.logCRB);
            acceptedThisInner = true;

            trustRadius = min(trustRadiusMax,trustExpand*trustRadius);
            rhoF = max(curvatureMin,rhoF/curvatureDecrease);
            rhoC = max(curvatureMin,rhoC/curvatureDecrease);
            break;
        end

        if acceptedThisInner
            break;
        end

        % Rebuild a more conservative model only after all cheap line-search
        % trials for the current direction have failed.
        rhoF = min(curvatureMax,curvatureGrowth*rhoF);
        rhoC = min(curvatureMax,curvatureGrowth*rhoC);
        trustRadius = max(trustRadiusMin,trustShrink*trustRadius);
    end

    if noUsefulDirection
        info.status = 'no useful feasible SCA phase direction.';
        break;
    end
    if ~acceptedThisInner
        info.status = ['phase retained; no exact feasible monotone step ', ...
            'was found within the model/line-search budget.'];
        break;
    end
end

[Pfinal,~] = build_P(theta,ch.Omega);
[Jfinal,crbFinalInfo] = SIM_crb_value(ch,Pfinal,Rx,alg);
[cFinal,~,~] = SIM_sinr_phase_constraints( ...
    params,ch,Pfinal,theta,Qcomm,Rs,phaseGammaTarget,constraintScale,false);

info.CRBEnd = Jfinal;
info.maxViolationEnd = max([0;cFinal]);
info.success = crbFinalInfo.valid && info.maxViolationEnd <= constraintTol && ...
    crbFinalInfo.logCRB <= crbInfoStart.logCRB+objectiveTol;
info.trustRadiusEnd = trustRadius;
info.curvatureFEnd = rhoF;
info.curvatureCMaxEnd = max(rhoC);
info.relMinEigAEnd = crbFinalInfo.relMinEigA;
info.relMinEigBEnd = crbFinalInfo.relMinEigB;
if strcmp(info.status,'not-started')
    if info.acceptedSteps > 0
        info.status = 'maximum phase SCA iterations reached.';
    else
        info.status = 'phase retained without an accepted step.';
    end
end
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
