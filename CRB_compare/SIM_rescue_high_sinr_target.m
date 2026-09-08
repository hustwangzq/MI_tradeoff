function [candidate,info] = SIM_rescue_high_sinr_target( ...
    params,ch,alg,baseInit,previousRecord,gammaValue,gammaTarget, ...
    commInfo,opts,sigma_c2,sigma_s2)
%SIM_RESCUE_HIGH_SINR_TARGET Strong rescue before terminating L=2/L=3.
%
% The rescue is invoked only after normal continuation/refinement fails.
% It tries the last valid Pareto state, the independently found pure-
% communication phase (if available), and small deterministic perturbations.
% Every trial still runs the original full CRB-SINR AO solver and is kept
% only when the exact RAW CRB and SINR post-checks succeed.

L = get_layer(ch);
candidate = make_candidate_template(gammaValue);
info = make_info_template(L,gammaValue);

if L < 2 || ~is_valid_record(previousRecord)
    info.status = 'high-SINR rescue not applicable';
    return;
end

maxOuter = select_by_layer(opts,'rescueMaxOuterByL',L,300);
maxPhase = select_by_layer(opts,'rescueMaxPhaseInnerByL',L,8);
if L >= 3
    maxOuter = max(maxOuter,360);
    maxPhase = max(maxPhase,10);
end
numPerturb = max(1,round(select_by_layer(opts,'rescueNumPerturbByL',L,2)));
perturbStd = select_by_layer(opts,'rescuePerturbStdByL',L,0.025);
seedBase = round(get_option(opts,'seed',930000)) + 10000*L + ...
    round(100*10*log10(max(gammaValue,realmin)));

algRescue = alg;
algRescue.maxOuter = maxOuter;
algRescue.maxPhaseInner = maxPhase;
algRescue.outerMinIterations = min(maxOuter,max(25, ...
    get_option(alg,'outerMinIterations',15)));
algRescue.outerStableIterations = max(8,get_option(alg,'outerStableIterations',6));
algRescue.useOuterStop = true;
algRescue.useRobustSDR = true;
algRescue.robustSDRAdaptive = true;
algRescue.robustSDRFixedFallbackRegRel = ...
    get_option(opts,'fixedFallbackRegRel',1e-7);
% Near the communication boundary allow raw matrices down to the existing
% physical hard threshold, but do not replace raw inverse checks.
algRescue.crbHardMinRelEig = min(get_option(alg,'crbHardMinRelEig',1e-12),1e-12);
algRescue.crbMinRcond = min(get_option(alg,'crbMinRcond',1e-12),1e-12);
if L >= 3
    algRescue.phaseTrustRadius0 = min(get_option(alg,'phaseTrustRadius0',0.08),0.03);
    algRescue.phaseTrustRadiusMax = min(get_option(alg,'phaseTrustRadiusMax',0.20),0.06);
else
    algRescue.phaseTrustRadius0 = min(get_option(alg,'phaseTrustRadius0',0.08),0.05);
    algRescue.phaseTrustRadiusMax = min(get_option(alg,'phaseTrustRadiusMax',0.20),0.10);
end
algRescue.verbose = false;
algRescue.progressEvery = max(maxOuter+1,1000000);

initList = cell(0,1);
labelList = cell(0,1);
basePrev = record_to_init(baseInit,previousRecord,sigma_c2,sigma_s2);
[initList,labelList] = append_unique_init(initList,labelList,basePrev,'previous-report');

% Communication-optimal phase provides a second branch seed, but the full
% CRB-SINR AO solver is rerun at the requested target; the communication
% solution itself is never written to the Pareto curve.
if isstruct(commInfo) && isfield(commInfo,'success') && commInfo.success && ...
        isfield(commInfo,'bestTheta') && ~isempty(commInfo.bestTheta)
    commInit = basePrev;
    commInit.theta0 = commInfo.bestTheta;
    [initList,labelList] = append_unique_init( ...
        initList,labelList,commInit,'pure-comm-phase');
end

oldStream = rng;
cleanupObj = onCleanup(@() rng(oldStream)); %#ok<NASGU>
rng(seedBase,'twister');
baseCount = numel(initList);
for ib = 1:baseCount
    for ip = 1:numPerturb
        pert = initList{ib};
        if isfield(pert,'theta0') && ~isempty(pert.theta0)
            pert.theta0 = pert.theta0 + perturbStd*randn(size(pert.theta0));
            [initList,labelList] = append_unique_init(initList,labelList,pert, ...
                sprintf('%s-pert-%d',labelList{ib},ip));
        end
    end
end

bestCRB = Inf;
for ii = 1:numel(initList)
    [solTry,histTry] = SIM_run_CRB_SINR_AO_solver( ...
        params,ch,algRescue,initList{ii},gammaTarget);
    info.numAttempts = info.numAttempts + 1;
    info.attemptLabels{end+1,1} = labelList{ii};
    info.attemptSuccess(end+1,1) = is_successful_solution(solTry);
    info.attemptStatus{end+1,1} = status_text(solTry,'status','solver-failed');

    if ~is_successful_solution(solTry)
        continue;
    end
    crbPhysical = (solTry.sigma_s2/params.T)*solTry.metrics.CRB;
    if ~(isfinite(crbPhysical) && crbPhysical > 0)
        continue;
    end
    if crbPhysical < bestCRB
        bestCRB = crbPhysical;
        candidate = make_candidate_from_solution( ...
            ['high-SINR-rescue-',labelList{ii}],solTry,histTry,gammaValue);
        info.bestLabel = labelList{ii};
        info.bestCRBPhysicalDB = 10*log10(crbPhysical);
        info.bestMinSINRdB = 10*log10(max(min(solTry.metrics.sinr),realmin));
    end
end

info.success = candidate.success;
if info.success
    info.status = 'high-SINR rescue found exact feasible complete solution';
else
    info.status = 'high-SINR rescue exhausted without exact feasible solution';
end
end

function candidate = make_candidate_template(gammaValue)
sol0 = struct('success',false,'status','not-run','metrics', ...
    struct('CRB',Inf,'sinr',zeros(0,1),'crbInfo', ...
    struct('relMinEigA',NaN,'relMinEigB',NaN)));
hist0 = struct('status','not-run');
candidate = struct( ...
    'label','none', ...
    'sol',sol0, ...
    'hist',hist0, ...
    'trialGamma',gammaValue, ...
    'selectable',true, ...
    'success',false, ...
    'CRB',Inf, ...
    'minSINR',NaN, ...
    'relMinEigA',NaN, ...
    'relMinEigB',NaN, ...
    'status','not-run');
end

function candidate = make_candidate_from_solution(label,sol,hist,gammaValue)
candidate = make_candidate_template(gammaValue);
candidate.label = label;
candidate.sol = sol;
candidate.hist = hist;
candidate.success = true;
candidate.CRB = sol.metrics.CRB;
candidate.minSINR = min(sol.metrics.sinr);
candidate.relMinEigA = sol.metrics.crbInfo.relMinEigA;
candidate.relMinEigB = sol.metrics.crbInfo.relMinEigB;
candidate.status = status_text(sol,'status','solved');
end

function info = make_info_template(L,gammaValue)
info = struct( ...
    'success',false, ...
    'L',L, ...
    'gammaValue',gammaValue, ...
    'gammaValueDB',10*log10(max(gammaValue,realmin)), ...
    'numAttempts',0, ...
    'attemptLabels',{cell(0,1)}, ...
    'attemptSuccess',false(0,1), ...
    'attemptStatus',{cell(0,1)}, ...
    'bestLabel','none', ...
    'bestCRBPhysicalDB',NaN, ...
    'bestMinSINRdB',NaN, ...
    'status','not-run');
end

function init = record_to_init(baseInit,record,sigma_c2,sigma_s2)
init = baseInit;
if isfield(record,'W') && ~isempty(record.W), init.W0 = record.W; end
if isfield(record,'QcommSDR') && ~isempty(record.QcommSDR), init.Qcomm0=record.QcommSDR; end
if isfield(record,'RsSDR') && ~isempty(record.RsSDR), init.Rs0=record.RsSDR; end
if isfield(record,'theta') && ~isempty(record.theta), init.theta0 = record.theta; end
if isfield(record,'S') && ~isempty(record.S), init.S = record.S; end
init.sigma_c2 = sigma_c2;
init.sigma_s2 = sigma_s2;
end

function [initList,labelList] = append_unique_init(initList,labelList,init,label)
if ~isstruct(init) || ~isfield(init,'theta0') || isempty(init.theta0)
    return;
end
for ii = 1:numel(initList)
    if isfield(initList{ii},'theta0') && ...
            isequal(size(initList{ii}.theta0),size(init.theta0))
        delta = angle(exp(1j*(initList{ii}.theta0-init.theta0)));
        if norm(delta(:)) <= 1e-10*sqrt(max(numel(delta),1))
            return;
        end
    end
end
initList{end+1,1} = init;
labelList{end+1,1} = label;
end

function tf = is_valid_record(record)
tf = isstruct(record) && isfield(record,'success') && ...
    isscalar(record.success) && logical(record.success) && ...
    isfield(record,'QcommSDR') && ~isempty(record.QcommSDR) && ...
    isfield(record,'theta') && ~isempty(record.theta);
end

function tf = is_successful_solution(sol)
tf = isstruct(sol) && isfield(sol,'success') && isscalar(sol.success) && ...
    logical(sol.success) && isfield(sol,'metrics') && isstruct(sol.metrics) && ...
    isfield(sol.metrics,'CRB') && isfinite(sol.metrics.CRB) && ...
    sol.metrics.CRB > 0 && ...
    logical(get_option(sol.metrics,'crbValid',false)) && ...
    logical(get_option(sol.metrics,'sinrFeasible',false));
end

function L = get_layer(ch)
L = 1;
if isfield(ch,'L') && isscalar(ch.L) && isfinite(ch.L)
    L = max(1,round(ch.L));
end
end

function value = select_by_layer(opts,name,L,defaultValue)
value = defaultValue;
if isstruct(opts) && isfield(opts,name) && ~isempty(opts.(name))
    x = opts.(name);
    if isnumeric(x) && numel(x) >= L && isfinite(x(L))
        value = x(L);
    end
end
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end

function textValue = status_text(s,fieldName,defaultValue)
textValue = defaultValue;
if isstruct(s) && isfield(s,fieldName) && ~isempty(s.(fieldName))
    textValue = char(string(s.(fieldName)));
end
end
