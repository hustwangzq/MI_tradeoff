function [carrier,info] = SIM_build_comm_carrier( ...
    params,chIn,alg,bestInit,gammaCurrent,gammaNext,previousCarrier)
%SIM_BUILD_COMM_CARRIER Build a hidden covariance-domain continuation state.
%
% bestInit remains the sensing/report state.  previousCarrier is independently
% propagated across formal SINR points.  The returned carrier is not itself
% plotted; it competes with bestInit only after a complete sensing AO solve.
% No rank-one recovery is performed anywhere in this routine.

if nargin<7, previousCarrier=[]; end

carrier = [];
info = struct('success',false,'status','not-started', ...
    'currentGamma',gammaCurrent,'nextGamma',gammaNext, ...
    'requiredCurrentMargin',gammaNext/max(gammaCurrent,realmin)-1, ...
    'startCurrentMargin',NaN,'outerCurrentMargin',NaN, ...
    'bestNextMargin',-Inf,'selectedMode','none','nextFeasible',false, ...
    'selectedCRBValid',false,'successReady',false,'reserveReady',false, ...
    'selectedClass',-Inf, ...
    'selectedCarrierCRB',Inf,'selectedCRBDeltaDB',Inf, ...
    'selectedQuality','none','selectedBlendWeight',NaN, ...
    'targetCurrentMargin',NaN,'targetNextMargin',NaN, ...
    'usedPreviousCarrier',false, ...
    'outerAcceptedSteps',0,'phaseIOuterAcceptedSteps',0, ...
    'phaseIAttempts',0,'allPhaseAcceptedSteps',0, ...
    'hiddenAcceptedSteps',0,'lastReachedGamma',gammaCurrent, ...
    'lastCVXStatus','not-run','goalFallbackUsed',false, ...
    'goalAttemptCount',0);

required = {'theta0','Qcomm0','Rs0','S','sigma_c2','sigma_s2'};
for ii = 1:numel(required)
    if ~isstruct(bestInit) || ~isfield(bestInit,required{ii}) || ...
            isempty(bestInit.(required{ii}))
        info.status = ['missing-',required{ii}];
        return;
    end
end
if ~(isscalar(gammaCurrent) && isfinite(gammaCurrent) && gammaCurrent>0 && ...
        isscalar(gammaNext) && isfinite(gammaNext) && gammaNext>gammaCurrent)
    info.status = 'invalid-gamma-transition';
    return;
end

requiredCurrentMargin=gammaNext/gammaCurrent-1;
if isfield(alg,'carrierCurrentTargetMarginRel') && ...
        ~isempty(alg.carrierCurrentTargetMarginRel)
    baseMargin=max(0,get_option(alg,'carrierCurrentTargetMarginRel',0.10));
    capMargin=max(baseMargin,get_option(alg, ...
        'carrierMaxCurrentTargetMarginRel',0.50));
    stepFactor=max(1,get_option(alg,'carrierStepReserveFactor',1.25));
    targetCurrentMargin=min(capMargin,max(baseMargin, ...
        stepFactor*requiredCurrentMargin));
    targetCurrentMargin=max(requiredCurrentMargin,targetCurrentMargin);
    gammaGoal=gammaCurrent*(1+targetCurrentMargin);
    deliveryMargin=max(0,gammaGoal/gammaNext-1);
else
    % Compatibility for previously saved local-validation settings.
    deliveryMargin=max(0,get_option(alg,'carrierDeliveryMarginRel',1e-2));
    gammaGoal=gammaNext*(1+deliveryMargin);
    targetCurrentMargin=gammaGoal/gammaCurrent-1;
end
info.targetCurrentMargin=targetCurrentMargin;
info.targetNextMargin=deliveryMargin;
reserveSatisfaction=max(0,min(1,get_option( ...
    alg,'carrierReserveSatisfactionRatio',0.85)));

ch = chIn;
ch.sigma_c2 = bestInit.sigma_c2;
ch.sigma_s2 = bestInit.sigma_s2;
gammaCurrentVector = gammaCurrent*ones(params.K,1);
gammaNextVector = gammaNext*ones(params.K,1);
theta = bestInit.theta0;
Qcomm = bestInit.Qcomm0;
Rs = bestInit.Rs0;
RxHealthReference = Rs;
for kk = 1:params.K
    RxHealthReference = RxHealthReference+Qcomm(:,:,kk);
end
[Pbest,~]=build_P(theta,ch.Omega);
[bestReferenceCRB,~]=SIM_crb_value(ch,Pbest,RxHealthReference,alg);

startMargin = min(relative_sinr_margin(theta,Qcomm,Rs,gammaCurrentVector));
info.startCurrentMargin = startMargin;
consider(theta,Qcomm,Rs,'bestCRB-seed');

% Stage A: Phi_L only.  This cannot change the extended-target CRB, and may
% already create enough reserve to cover the actual next target.
algOuter = alg;
algOuter.outerMarginMaxIter = max(round(get_option(alg, ...
    'carrierOuterMaxIter',30)),get_option(alg,'outerMarginMaxIter',12));
[thetaOuter,outerInfo] = SIM_update_outer_phase_margin( ...
    params,ch,theta,Qcomm,Rs,gammaCurrentVector,algOuter);
info.outerAcceptedSteps = get_option(outerInfo,'acceptedSteps',0);
outerMargin = min(relative_sinr_margin( ...
    thetaOuter,Qcomm,Rs,gammaCurrentVector));
info.outerCurrentMargin = outerMargin;
consider(thetaOuter,Qcomm,Rs,'outer-only');

priorSeedAvailable=has_complete_init(previousCarrier,params,size(theta))&& ...
    ~same_carrier_state(previousCarrier,theta,Qcomm,Rs);
thetaPriorOuter=[];
if priorSeedAvailable
    info.usedPreviousCarrier=true;
    consider(previousCarrier.theta0,previousCarrier.Qcomm0, ...
        previousCarrier.Rs0,'propagated-carrier');
    [thetaPriorOuter,priorOuterInfo]=SIM_update_outer_phase_margin( ...
        params,ch,previousCarrier.theta0,previousCarrier.Qcomm0, ...
        previousCarrier.Rs0,gammaCurrentVector,algOuter);
    info.outerAcceptedSteps=info.outerAcceptedSteps+ ...
        get_option(priorOuterInfo,'acceptedSteps',0);
    consider(thetaPriorOuter,previousCarrier.Qcomm0, ...
        previousCarrier.Rs0,'propagated-carrier-outer');
end

if info.reserveReady
    info.success = true;
    info.status = 'outer-only-next-feasible';
    return;
end

% Stage B: communication-only hidden-Gamma continuation.  Try the actual
% next target first.  If Phase-I cannot reach it, halve the Gamma step and
% let every feasible hidden point become the next carrier anchor.  CRB is
% deliberately not an acceptance requirement for hidden states.
maxPhaseI = max(1,round(get_option(alg,'carrierPhaseIMaxIter',32)));
minStepDB = get_option(alg,'carrierGammaMinStepDB',1e-3);
feasTol = get_option(alg,'carrierSINRRelativeTol',2e-5);
numRequestedSeeds = max(1,round(get_option(alg,'carrierNumPhaseSeeds',1)));
thetaSeeds = {thetaOuter};
seedLabels = {'best-outer'};
seedQ={Qcomm};
seedRs={Rs};
if priorSeedAvailable
    thetaSeeds{end+1}=thetaPriorOuter;
    seedLabels{end+1}='propagated-carrier';
    seedQ{end+1}=previousCarrier.Qcomm0;
    seedRs{end+1}=previousCarrier.Rs0;
end
if numRequestedSeeds > 1
    oldRng = rng;
    rng(round(get_option(alg,'carrierSeed',24681357))+1000*ch.L,'twister');
    perturbScale = get_option(alg,'carrierPhasePerturbScales',[0.05,0.20]);
    perturbScale = perturbScale(:).';
    for is = 2:numRequestedSeeds
        if is-1 <= numel(perturbScale)
            thetaSeeds{end+1} = mod(thetaOuter+ ...
                perturbScale(is-1)*randn(size(thetaOuter)),2*pi);
            seedLabels{end+1} = sprintf('perturb-%.3g',perturbScale(is-1));
        else
            thetaSeeds{end+1} = 2*pi*rand(size(thetaOuter));
            seedLabels{end+1} = sprintf('random-%d',is-numel(perturbScale)-1);
        end
        seedQ{end+1}=Qcomm;
        seedRs{end+1}=Rs;
    end
    rng(oldRng);
end

fallbackMargin=max(0,get_option(alg,'carrierFallbackDeliveryMarginRel',1e-2));
goalList=unique([gammaGoal,gammaNext*(1+fallbackMargin),gammaNext],'stable');
for ig=1:numel(goalList)
    activeGoal=goalList(ig);
    info.goalAttemptCount=ig;
    for is = 1:numel(thetaSeeds)
        thetaWork = thetaSeeds{is};
        Qwork = seedQ{is};
        RsWork = seedRs{is};
        gammaReached = gammaCurrent;
        gammaTrial = activeGoal;
        for it = 1:maxPhaseI
            info.phaseIAttempts = info.phaseIAttempts+1;
            gammaTrialVector = gammaTrial*ones(params.K,1);
            marginSol = SIM_solve_comm_margin_SDR( ...
                params,ch,thetaWork,gammaTrialVector,alg);
            info.lastCVXStatus = marginSol.cvxStatus;
            if ~marginSol.success, break; end
            Qwork = marginSol.Qcomm;
            RsWork = marginSol.Rs;
            consider(thetaWork,Qwork,RsWork, ...
                ['phase-I-covariance-',seedLabels{is}]);
            trialMargin = min(relative_sinr_margin( ...
                thetaWork,Qwork,RsWork,gammaTrialVector));

        % Re-optimize Phi_L after every covariance update.  This is the
        % central carrier alternation: Phi_L can improve communication
        % reserve without changing the extended-target CRB geometry.
            [thetaWork,outerPhaseInfo] = SIM_update_outer_phase_margin( ...
                params,ch,thetaWork,Qwork,RsWork,gammaTrialVector,algOuter);
            info.phaseIOuterAcceptedSteps = info.phaseIOuterAcceptedSteps+ ...
                get_option(outerPhaseInfo,'acceptedSteps',0);
            consider(thetaWork,Qwork,RsWork, ...
                ['phase-I-outer-',seedLabels{is}]);
            trialMargin = min(relative_sinr_margin( ...
                thetaWork,Qwork,RsWork,gammaTrialVector));

            if trialMargin < -feasTol
                algPhase = alg;
                algPhase.carrierCRBHealthRx = RxHealthReference;
                [thetaWork,phaseInfo] = SIM_update_all_phase_margin( ...
                    params,ch,thetaWork,Qwork,RsWork,gammaTrialVector,algPhase);
                info.allPhaseAcceptedSteps = info.allPhaseAcceptedSteps+ ...
                    get_option(phaseInfo,'acceptedSteps',0);
                consider(thetaWork,Qwork,RsWork, ...
                    ['phase-I-all-phase-',seedLabels{is}]);
                trialMargin = min(relative_sinr_margin( ...
                    thetaWork,Qwork,RsWork,gammaTrialVector));
            end

            if trialMargin >= -feasTol
                gammaReached = gammaTrial;
                info.hiddenAcceptedSteps = info.hiddenAcceptedSteps+1;
                info.lastReachedGamma = max(info.lastReachedGamma,gammaReached);
                if gammaReached >= activeGoal*(1-1e-12), break; end
                gammaTrial = activeGoal;
                continue;
            end

            stepDB = 10*log10(gammaTrial/max(gammaReached,realmin));
            if stepDB <= minStepDB*(1+1e-8), break; end
            gammaTrial = sqrt(gammaReached*gammaTrial);
        end
        if info.reserveReady, break; end
    end
    if info.selectedClass>=2, break; end
end
info.goalFallbackUsed=info.goalAttemptCount>1;

if isempty(carrier)
    info.status = 'no-usable-carrier';
elseif info.reserveReady
    info.success = true;
    info.successReady = true;
    info.status = 'phase-I-next-feasible';
elseif info.selectedClass >= 2
    % A strictly feasible full covariance/phase pair is still a valid
    % carrier even when the requested extra reserve was not fully reached.
    info.success = true;
    info.successReady = true;
    if info.goalFallbackUsed
        info.status='phase-I-next-feasible-reduced-reserve';
    else
        info.status='phase-I-next-feasible-partial-reserve';
    end
else
    % A partial carrier is useful only diagnostically; do not advertise it as
    % a valid next-target warm start.
    info.status = 'best-carrier-below-next-target';
end

    function consider(thetaCandidate,Qcandidate,Rscandidate,mode)
        [Pcandidate,~] = build_P(thetaCandidate,ch.Omega);
        blendList = unique([0,get_option(alg,'carrierBlendWeights', ...
            logspace(-6,0,25))]);
        blendList = blendList(isfinite(blendList)&blendList>=0&blendList<=1);
        formalMarginFloor = get_option(alg,'carrierFormalMarginFloor',0);
        for iblend = 1:numel(blendList)
            lambda = blendList(iblend);
            Qtrial = (1-lambda)*Qcandidate+lambda*Qcomm;
            Rstrial = (1-lambda)*Rscandidate+lambda*Rs;
            nextMargin = min(relative_sinr_margin( ...
                thetaCandidate,Qtrial,Rstrial,gammaNextVector));
            if ~isfinite(nextMargin), continue; end
            nextFeasible = nextMargin >= formalMarginFloor;
            % Avoid rebuilding an already useful propagated carrier merely
            % because it misses the nominal reserve by a small fraction.
            reserveFeasible = nextMargin >= ...
                reserveSatisfaction*deliveryMargin;
            RxTrial = Rstrial;
            for jj=1:params.K, RxTrial=RxTrial+Qtrial(:,:,jj); end
            [crbValue,crbCandidate] = SIM_crb_value( ...
                ch,Pcandidate,RxTrial,alg);
            crbValid = logical(get_option(crbCandidate,'valid',false));
            crbDeltaDB=Inf;
            if crbValid&&isfinite(bestReferenceCRB)&&bestReferenceCRB>0
                crbDeltaDB=10*log10(crbValue/bestReferenceCRB);
            end
            protectedLoss=get_option(alg,'carrierMaxCRBLossDB',6);
            relaxedLoss=max(protectedLoss,get_option( ...
                alg,'carrierRelaxedMaxCRBLossDB',15));
            quality='fallback';
            if nextFeasible&&crbValid
                if crbDeltaDB<=protectedLoss
                    quality='protected';
                    candidateClass=5+double(reserveFeasible);
                elseif crbDeltaDB<=relaxedLoss
                    quality='relaxed';
                    candidateClass=3+double(reserveFeasible);
                else
                    candidateClass=2;
                end
            else
                candidateClass=double(nextFeasible);
            end
            betterClass = candidateClass > info.selectedClass;
            sameClass=info.selectedClass==candidateClass;
            materiallyMoreMargin=nextMargin>info.bestNextMargin+ ...
                max(1e-6,0.02*abs(info.bestNextMargin));
            betterProtectedMargin=sameClass&&candidateClass>=3&& ...
                materiallyMoreMargin;
            comparableMargin=abs(nextMargin-info.bestNextMargin)<= ...
                max(1e-6,0.02*abs(info.bestNextMargin));
            betterCRB=sameClass&&candidateClass>=2&& ...
                (candidateClass==2||comparableMargin)&& ...
                crbValue<info.selectedCarrierCRB*(1-1e-10);
            betterMargin = candidateClass == info.selectedClass && ...
                candidateClass < 2 && nextMargin > info.bestNextMargin+1e-12;
            if ~(betterClass||betterProtectedMargin||betterCRB||betterMargin)
                continue;
            end
            candidate = bestInit;
            candidate.theta0 = thetaCandidate;
            candidate.Qcomm0 = Qtrial;
            candidate.Rs0 = Rstrial;
            candidate.carrierMode = char(mode);
            candidate.carrierNextMargin = nextMargin;
            candidate.carrierTargetGamma = gammaNext;
            candidate.carrierBlendWeight = lambda;
            candidate.carrierCRBDeltaDB=crbDeltaDB;
            carrier = candidate;
            info.bestNextMargin = nextMargin;
            info.selectedMode = char(mode);
            info.nextFeasible = nextFeasible;
            info.selectedCRBValid = crbValid;
            info.reserveReady = reserveFeasible&&crbValid;
            info.successReady = info.reserveReady;
            info.selectedClass = candidateClass;
            info.selectedCarrierCRB = crbValue;
            info.selectedCRBDeltaDB=crbDeltaDB;
            info.selectedQuality=quality;
            info.selectedBlendWeight = lambda;
        end
    end

    function relativeMargin = relative_sinr_margin( ...
            thetaValue,Qvalue,RsValue,gammaValue)
        [Pvalue,~] = build_P(thetaValue,ch.Omega);
        RxValue = RsValue;
        for jj = 1:params.K, RxValue=RxValue+Qvalue(:,:,jj); end
        relativeMargin = zeros(params.K,1);
        for kk = 1:params.K
            gk = ch.G'*Pvalue'*ch.hUsers{kk};
            desired = real(gk'*Qvalue(:,:,kk)*gk);
            denominator = real(gk'*(RxValue-Qvalue(:,:,kk))*gk)+ch.sigma_c2;
            sinrValue = desired/max(denominator,realmin);
            relativeMargin(kk) = sinrValue/max(gammaValue(kk),realmin)-1;
        end
    end
end

function tf=has_complete_init(init,params,thetaSize)
tf=isstruct(init)&&isfield(init,'theta0')&& ...
    isequal(size(init.theta0),thetaSize)&& ...
    isfield(init,'Qcomm0')&& ...
    isequal(size(init.Qcomm0),[params.Nt,params.Nt,params.K])&& ...
    isfield(init,'Rs0')&&isequal(size(init.Rs0),[params.Nt,params.Nt]);
end

function tf=same_carrier_state(init,theta,Qcomm,Rs)
phaseSame=norm(exp(1j*init.theta0(:))-exp(1j*theta(:)))<= ...
    1e-12*sqrt(numel(theta));
qScale=max([norm(init.Qcomm0(:)),norm(Qcomm(:)),realmin]);
rScale=max([norm(init.Rs0(:)),norm(Rs(:)),realmin]);
tf=phaseSame&&norm(init.Qcomm0(:)-Qcomm(:))<=1e-10*qScale&& ...
    norm(init.Rs0(:)-Rs(:))<=1e-10*rScale;
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
