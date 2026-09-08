function sdr = SIM_solve_CRB_SINR_SDR_robust(params,ch,P,gammaTarget,alg)
%SIM_SOLVE_CRB_SINR_SDR_ROBUST Adaptive robust wrapper for SIM L=2/3.
%
% Adaptive mode tries increasingly strong INTERNAL A/B regularization while
% relaxing the optional sensing-eigenvalue floor. The first candidate that
% passes exact raw SINR, PSD, power, and unregularized CRB checks is returned.
% If adaptive loading does not help, the final attempt uses a configurable
% fixed regularization (default 1e-7), as requested for a conservative
% fallback.

if nargin < 5 || isempty(alg)
    alg = struct('useRobustSDR',true);
end

useRobust = logical(get_option(alg,'useRobustSDR',true));
% Keep L=1 on the historical single-SDR path. Robust/adaptive retries are
% reserved for the numerically harder L=2/L=3 cases.
layerValue = 1;
if isstruct(ch) && isfield(ch,'L') && isscalar(ch.L) && isfinite(ch.L)
    layerValue = max(1,round(ch.L));
end
if ~useRobust || layerValue <= 1
    sdr = SIM_solve_CRB_SINR_SDR(params,ch,P,gammaTarget,alg);
    sdr = attach_robust_fields_single(sdr,alg);
    return;
end

adaptiveMode = logical(get_option(alg,'robustSDRAdaptive',true));
if adaptiveMode
    plan = SIM_build_adaptive_CRB_regularization_plan(ch,P,alg);
else
    plan = build_legacy_plan(P,alg);
end

numAttempts = plan.numAttempts;
attemptStatus = repmat({'not-run'},1,numAttempts);
attemptSlvtol = NaN(1,numAttempts);
attemptRegA = NaN(1,numAttempts);
attemptRegB = NaN(1,numAttempts);
attemptFloor = NaN(1,numAttempts);
attemptScale = NaN(1,numAttempts);
attemptPrecision = repmat({'default'},1,numAttempts);
selected = 0;

% First attempt initializes a fully populated base-SDR schema.
algTry = make_attempt_alg(alg,plan,1);
sdr = SIM_solve_CRB_SINR_SDR(params,ch,P,gammaTarget,algTry);

for ia = 1:numAttempts
    if ia > 1
        algTry = make_attempt_alg(alg,plan,ia);
        try
            cvx_clear;
        catch
        end
        sdrTry = SIM_solve_CRB_SINR_SDR(params,ch,P,gammaTarget,algTry);
    else
        sdrTry = sdr;
    end

    attemptStatus{ia} = sdrTry.status;
    attemptSlvtol(ia) = safe_scalar(sdrTry,'cvxSlvtol',NaN);
    attemptRegA(ia) = plan.regRelA(ia);
    attemptRegB(ia) = plan.regRelB(ia);
    attemptFloor(ia) = plan.sensingFloorRel(ia);
    attemptScale(ia) = plan.schurScale(ia);
    attemptPrecision{ia} = plan.precision{ia};
    sdr = sdrTry;

    if sdrTry.success
        selected = ia;
        break;
    end
end

attemptsUsed = find_last_attempt(attemptStatus);
if selected == 0
    selected = attemptsUsed;
    sdr.status = sprintf(['Adaptive robust SDR exhausted %d attempt(s). ', ...
        'Last: %s'],attemptsUsed,sdr.status);
end

sdr.robustUsed = true;
sdr.robustAdaptive = adaptiveMode;
sdr.robustPlanMode = plan.mode;
sdr.robustAttemptCount = attemptsUsed;
sdr.robustSelectedAttempt = selected;
sdr.robustSelectedRegRelA = attemptRegA(selected);
sdr.robustSelectedRegRelB = attemptRegB(selected);
sdr.robustSelectedRegRel = max(attemptRegA(selected),attemptRegB(selected));
sdr.robustSelectedSensingFloorRel = attemptFloor(selected);
sdr.robustSelectedPrecision = attemptPrecision{selected};
sdr.robustSelectedScale = attemptScale(selected);
sdr.robustAttemptStatus = attemptStatus(1:attemptsUsed);
sdr.robustAttemptSlvtol = attemptSlvtol(1:attemptsUsed);
sdr.robustAttemptRegRelA = attemptRegA(1:attemptsUsed);
sdr.robustAttemptRegRelB = attemptRegB(1:attemptsUsed);
sdr.robustAttemptRegRel = max(attemptRegA(1:attemptsUsed), ...
    attemptRegB(1:attemptsUsed));
sdr.robustAttemptSensingFloorRel = attemptFloor(1:attemptsUsed);
sdr.robustAttemptScale = attemptScale(1:attemptsUsed);
sdr.robustAttemptPrecision = attemptPrecision(1:attemptsUsed);
sdr.robustRawRelMinEigBAtPlan = plan.rawRelMinEigB;
end

function algTry = make_attempt_alg(alg,plan,ia)
algTry = alg;
algTry.cvxSchurScale = plan.schurScale(ia);
algTry.cvxPrecisionMode = char(string(plan.precision{ia}));
algTry.useInternalCRBRegularization = ...
    plan.regRelA(ia) > 0 || plan.regRelB(ia) > 0;
algTry.crbInternalRegRelA = plan.regRelA(ia);
algTry.crbInternalRegRelB = plan.regRelB(ia);
% Compatibility alias for older phase-gradient code.
algTry.crbInternalRegRel = max(plan.regRelA(ia),plan.regRelB(ia));
algTry.sensingEigFloorRel = plan.sensingFloorRel(ia);
end

function plan = build_legacy_plan(P,alg)
N = size(P,1);
scaleValue = get_option(alg,'robustSDRSchurScale',N);
if ischar(scaleValue) || isstring(scaleValue)
    if strcmpi(char(string(scaleValue)),'N')
        scaleValue = N;
    else
        scaleValue = N;
    end
end
if ~(isnumeric(scaleValue) && isscalar(scaleValue) && ...
        isfinite(scaleValue) && scaleValue > 0)
    scaleValue = N;
end
regList = get_option(alg,'robustSDRRegList',[0,1e-8,1e-7]);
regList = regList(:).';
regList = regList(isfinite(regList) & regList >= 0);
if isempty(regList)
    regList = [0,1e-8,1e-7];
end
precision = repmat({'default'},1,numel(regList));
plan = struct( ...
    'numAttempts',numel(regList), ...
    'regRelA',regList, ...
    'regRelB',regList, ...
    'sensingFloorRel',get_option(alg,'sensingEigFloorRel',0)*ones(size(regList)), ...
    'precision',{precision}, ...
    'schurScale',scaleValue*ones(size(regList)), ...
    'rawRelMinEigB',NaN, ...
    'fixedFallbackRegRel',regList(end), ...
    'mode','legacy-fixed-list');
end

function sdr = attach_robust_fields_single(sdr,alg)
regA = get_option(alg,'crbInternalRegRelA', ...
    get_option(alg,'crbInternalRegRel',0));
regB = get_option(alg,'crbInternalRegRelB', ...
    get_option(alg,'crbInternalRegRel',0));
sdr.robustUsed = false;
sdr.robustAdaptive = false;
sdr.robustPlanMode = 'single-base-solve';
sdr.robustAttemptCount = 1;
sdr.robustSelectedAttempt = 1;
sdr.robustSelectedRegRelA = regA;
sdr.robustSelectedRegRelB = regB;
sdr.robustSelectedRegRel = max(regA,regB);
sdr.robustSelectedSensingFloorRel = get_option(alg,'sensingEigFloorRel',0);
sdr.robustSelectedPrecision = char(string(get_option(alg,'cvxPrecisionMode','high')));
sdr.robustSelectedScale = get_option(alg,'cvxSchurScale',1);
sdr.robustAttemptStatus = {sdr.status};
sdr.robustAttemptSlvtol = safe_scalar(sdr,'cvxSlvtol',NaN);
sdr.robustAttemptRegRelA = regA;
sdr.robustAttemptRegRelB = regB;
sdr.robustAttemptRegRel = max(regA,regB);
sdr.robustAttemptSensingFloorRel = sdr.robustSelectedSensingFloorRel;
sdr.robustAttemptScale = sdr.robustSelectedScale;
sdr.robustAttemptPrecision = {sdr.robustSelectedPrecision};
sdr.robustRawRelMinEigBAtPlan = safe_scalar(sdr,'relMinEigB',NaN);
end

function idx = find_last_attempt(statusList)
idx = 1;
for ii = 1:numel(statusList)
    if ~strcmp(statusList{ii},'not-run')
        idx = ii;
    end
end
end

function value = safe_scalar(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && isscalar(s.(name)) && ...
        isnumeric(s.(name))
    value = s.(name);
end
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
