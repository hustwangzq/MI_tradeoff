function [gammaCommMax,info] = SIM_estimate_pure_comm_sinr_boundary( ...
    params,ch,initList,opts)
%SIM_ESTIMATE_PURE_COMM_SINR_BOUNDARY Estimate pure-communication ceiling.
%
% If opts.referenceDBByL(L) is finite, that value is used directly. This is
% the preferred mode when an independent pure-communication simulation has
% already established the corresponding limit (for example ~29 dB).
%
% Otherwise, a reproducible multi-start fixed-phase communication-only SDR
% search is performed. Its result is a numerical lower bound/estimate, not a
% proof of the global joint phase/beamforming optimum.

if nargin < 4 || isempty(opts)
    opts = struct('enable',true);
end

L = get_layer(ch);
manualDB = NaN;
if isfield(opts,'referenceDBByL') && numel(opts.referenceDBByL) >= L
    manualDB = opts.referenceDBByL(L);
end
if isfinite(manualDB)
    gammaCommMax = 10^(manualDB/10);
    info = make_info_template(L);
    info.success = true;
    info.gammaCommMax = gammaCommMax;
    info.gammaCommMaxDB = manualDB;
    info.source = 'manual-pure-communication-reference';
    info.status = 'manual reference used';
    return;
end

if ~logical(get_option(opts,'enable',true))
    gammaCommMax = NaN;
    info = make_info_template(L);
    info.status = 'communication-boundary-estimator-disabled';
    return;
end

seedBase = round(get_option(opts,'seed',910000)) + 1000*L;
numRandom = max(0,round(select_by_layer(opts,'numRandomByL',L,6)));
maxCandidateCap = max(1,round(select_by_layer(opts,'maxCandidatesByL',L,12)));

% Candidate phases are stored as a cell array to avoid structure-array
% schema issues. Existing initialization phases are tried first.
thetaList = cell(0,1);
labelList = cell(0,1);
for ii = 1:numel(initList)
    init = initList{ii};
    if isstruct(init) && isfield(init,'theta0') && ~isempty(init.theta0)
        [thetaList,labelList] = append_unique_theta(thetaList,labelList, ...
            init.theta0,sprintf('existing-init-%02d',ii));
    end
    if numel(thetaList) >= maxCandidateCap
        break;
    end
end

% Add deterministic random phases if needed.
oldStream = rng;
cleanupObj = onCleanup(@() rng(oldStream)); %#ok<NASGU>
rng(seedBase,'twister');
for ir = 1:numRandom
    if numel(thetaList) >= maxCandidateCap
        break;
    end
    thetaRand = 2*pi*rand(params.N,L);
    [thetaList,labelList] = append_unique_theta(thetaList,labelList, ...
        thetaRand,sprintf('random-%02d',ir));
end

if isempty(thetaList)
    thetaList{1,1} = zeros(params.N,L);
    labelList{1,1} = 'zero-phase-fallback';
end

solverOpts = struct( ...
    'gammaSeed',get_option(opts,'gammaSeed',1), ...
    'maxExpand',get_option(opts,'maxExpand',35), ...
    'maxBisection',get_option(opts,'maxBisection',32), ...
    'upperCapDB',get_option(opts,'upperCapDB',45), ...
    'relativeTolerance',get_option(opts,'relativeTolerance',2e-4), ...
    'cvxQuiet',get_option(opts,'cvxQuiet',true), ...
    'cvxPrecisionMode',get_option(opts,'cvxPrecisionMode','default'));
if isfield(opts,'cvxSolver') && ~isempty(opts.cvxSolver)
    solverOpts.cvxSolver = opts.cvxSolver;
end

gammaCommMax = NaN;
bestTheta = thetaList{1};
bestLabel = 'none';
bestSol = make_empty_comm_sol(params);
trialGammaDB = NaN(numel(thetaList),1);
trialSuccess = false(numel(thetaList),1);
trialStatus = repmat({''},numel(thetaList),1);

for ii = 1:numel(thetaList)
    [P,~] = build_P(thetaList{ii},ch.Omega);
    [gammaTmp,solTmp,infoTmp] = SIM_max_common_sinr_fixed_phase_SDR( ...
        params,ch,P,solverOpts);
    trialGammaDB(ii) = 10*log10(max(gammaTmp,realmin));
    trialSuccess(ii) = isfinite(gammaTmp) && infoTmp.success;
    trialStatus{ii} = infoTmp.status;
    fprintf('  [pure-comm L=%d %02d/%02d] %-20s  maxSINR %.3f dB\n', ...
        L,ii,numel(thetaList),labelList{ii},trialGammaDB(ii));
    if trialSuccess(ii) && (~isfinite(gammaCommMax) || gammaTmp > gammaCommMax)
        gammaCommMax = gammaTmp;
        bestTheta = thetaList{ii};
        bestLabel = labelList{ii};
        bestSol = solTmp;
    end
end

info = make_info_template(L);
info.success = isfinite(gammaCommMax) && gammaCommMax > 0;
info.gammaCommMax = gammaCommMax;
info.gammaCommMaxDB = 10*log10(max(gammaCommMax,realmin));
info.bestTheta = bestTheta;
info.bestLabel = bestLabel;
info.bestCommSolution = bestSol;
info.numCandidates = numel(thetaList);
info.trialLabels = labelList;
info.trialGammaMaxDB = trialGammaDB;
info.trialSuccess = trialSuccess;
info.trialStatus = trialStatus;
info.source = 'multistart-fixed-phase-pure-communication-SDR';
if info.success
    info.status = 'automatic pure-communication estimate available';
else
    info.status = 'automatic pure-communication estimate failed';
end
end

function info = make_info_template(L)
emptySol = struct('success',false,'gammaMax',NaN,'gammaMaxDB',NaN, ...
    'Qcomm',zeros(0,0,0),'sinr',zeros(0,1),'status','not-run');
info = struct( ...
    'success',false, ...
    'L',L, ...
    'gammaCommMax',NaN, ...
    'gammaCommMaxDB',NaN, ...
    'bestTheta',zeros(0,0), ...
    'bestLabel','none', ...
    'bestCommSolution',emptySol, ...
    'numCandidates',0, ...
    'trialLabels',{cell(0,1)}, ...
    'trialGammaMaxDB',zeros(0,1), ...
    'trialSuccess',false(0,1), ...
    'trialStatus',{cell(0,1)}, ...
    'source','none', ...
    'status','not-run');
end

function sol = make_empty_comm_sol(params)
sol = struct('success',false,'gammaMax',NaN,'gammaMaxDB',NaN, ...
    'Qcomm',zeros(params.Nt,params.Nt,params.K), ...
    'sinr',zeros(params.K,1),'status','not-run');
end

function [thetaList,labelList] = append_unique_theta(thetaList,labelList,theta,label)
if isempty(theta) || any(~isfinite(theta(:)))
    return;
end
for ii = 1:numel(thetaList)
    if isequal(size(thetaList{ii}),size(theta))
        delta = angle(exp(1j*(thetaList{ii}-theta)));
        if norm(delta(:)) <= 1e-10*sqrt(max(numel(delta),1))
            return;
        end
    end
end
thetaList{end+1,1} = theta;
labelList{end+1,1} = label;
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
