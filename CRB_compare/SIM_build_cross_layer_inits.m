function [initList,labels] = SIM_build_cross_layer_inits( ...
    params,chCurrent,baseInit,previousSol,opts)
%SIM_BUILD_CROSS_LAYER_INITS SIM initialization builder.
%
% Profiles:
%   compact    : at most two candidates, intended for anchor/recovery use.
%   multistart : deterministic base/structured candidates followed by
%                independent random-phase candidates up to maxCandidates.
%
% Cross-layer inheritance transfers phase only. W is never inherited across
% different L values because the current physical channel is different.

if nargin < 5 || isempty(opts)
    opts = struct();
end
opts = fill_defaults(opts,params);

L = chCurrent.L;
N = params.N;
if ~ismember(L,[1,2,3])
    error('SIM_build_cross_layer_inits:UnsupportedLayerCount', ...
        'The builder currently supports L=1,2,3.');
end
if ~isfield(baseInit,'theta0') || ~isequal(size(baseInit.theta0),[N,L])
    error('SIM_build_cross_layer_inits:InvalidBaseTheta', ...
        'baseInit.theta0 must have size N-by-chCurrent.L.');
end
if ~any(strcmpi(opts.profile,{'compact','multistart'}))
    error('SIM_build_cross_layer_inits:UnsupportedProfile', ...
        'opts.profile must be ''compact'' or ''multistart''.');
end

rngState = rng;
cleanupObj = onCleanup(@()rng(rngState)); %#ok<NASGU>
rng(opts.seed,'twister');

initList = {};
labels = {};

[initStructured,structuredLabel] = make_structured_start( ...
    baseInit,previousSol,N,L);

if strcmpi(opts.profile,'compact')
    % Candidate 1: random-phase warm start.
    initRandom = make_random_start(params,chCurrent,baseInit,opts,1,N,L);
    [initList,labels] = append_candidate(initList,labels,initRandom, ...
        'current-layer-random-phase-warm-start');

    % Candidate 2: deterministic zero/cross-layer structured start.
    [initList,labels] = append_candidate(initList,labels,initStructured, ...
        structuredLabel);
else
    % Multistart range search: first test the actual current-layer base phase,
    % then the structured cross-layer/zero phase, and finally independent
    % random phases. This keeps range-search diversity out of the formal
    % CRB continuation path.
    initBaseCandidate = copy_common_fields(baseInit,baseInit);
    [initList,labels] = append_candidate(initList,labels, ...
        initBaseCandidate,'current-layer-base-phase');
    [initList,labels] = append_candidate(initList,labels, ...
        initStructured,structuredLabel);

    randomIndex = 1;
    while numel(initList) < opts.maxCandidates
        initRandom = make_random_start(params,chCurrent,baseInit,opts, ...
            randomIndex,N,L);
        [initList,labels] = append_candidate(initList,labels,initRandom, ...
            sprintf('range-random-phase-%02d',randomIndex));
        randomIndex = randomIndex + 1;
        if randomIndex > 10*max(1,opts.maxCandidates)
            break;
        end
    end
end

[initList,labels] = limit_candidates(initList,labels,opts.maxCandidates);
end

function initRandom = make_random_start( ...
    params,chCurrent,baseInit,opts,randomIndex,N,L)
if opts.randomizeW
    initRandom = initialize_solver_state( ...
        params,chCurrent,opts.seed+100+randomIndex);
    initRandom = copy_stream_and_noise(initRandom,baseInit);
else
    initRandom = baseInit;
end
initRandom.theta0 = 2*pi*rand(N,L);
initRandom = copy_common_fields(initRandom,baseInit);
end

function [initStructured,label] = make_structured_start( ...
    baseInit,previousSol,N,L)
[previousTheta,previousValid] = extract_previous_theta(previousSol,N,L-1);
initStructured = baseInit;
if L == 1 || ~previousValid
    initStructured.theta0 = zeros(N,L);
    if L == 1
        label = 'current-layer-zero-phase';
    else
        label = 'current-layer-zero-phase-no-valid-previous-layer';
    end
else
    thetaStructured = zeros(N,L);
    if L == 2
        thetaStructured(:,2) = previousTheta(:,1);
        label = 'structured-inherit-L1-outer-to-L2-outer';
    else
        thetaStructured(:,2) = previousTheta(:,1);
        thetaStructured(:,3) = previousTheta(:,2);
        label = 'structured-inherit-L2-shift-outward-to-L3';
    end
    initStructured.theta0 = mod(thetaStructured,2*pi);
end
initStructured = copy_common_fields(initStructured,baseInit);
end

function [theta,valid] = extract_previous_theta(previousSol,N,expectedL)
theta = [];
valid = false;
if expectedL < 1 || isempty(previousSol) || ~isstruct(previousSol) || ...
        ~isfield(previousSol,'theta') || isempty(previousSol.theta)
    return;
end
theta = previousSol.theta;
valid = isequal(size(theta),[N,expectedL]) && all(isfinite(theta(:)));
if ~valid
    theta = [];
end
end

function initOut = copy_stream_and_noise(initOut,baseInit)
fields = {'S','sigma_c2','sigma_s2'};
for ii = 1:numel(fields)
    name = fields{ii};
    if isfield(baseInit,name) && ~isempty(baseInit.(name))
        initOut.(name) = baseInit.(name);
    end
end
end

function initOut = copy_common_fields(initOut,baseInit)
fields = {'W0','S','sigma_c2','sigma_s2'};
for ii = 1:numel(fields)
    name = fields{ii};
    if isfield(baseInit,name) && ~isempty(baseInit.(name))
        initOut.(name) = baseInit.(name);
    end
end
end

function [list,labels] = append_candidate(list,labels,candidate,label)
for ii = 1:numel(list)
    if isequal(size(list{ii}.theta0),size(candidate.theta0)) && ...
            norm(exp(1j*list{ii}.theta0(:))-exp(1j*candidate.theta0(:))) ...
            <= 1e-12*sqrt(numel(candidate.theta0))
        return;
    end
end
list{end+1,1} = candidate; %#ok<AGROW>
labels{end+1,1} = char(string(label)); %#ok<AGROW>
end

function [list,labels] = limit_candidates(list,labels,maxCandidates)
if isfinite(maxCandidates) && numel(list) > maxCandidates
    list = list(1:maxCandidates);
    labels = labels(1:maxCandidates);
end
end

function opts = fill_defaults(opts,params)
defaults.seed = params.seed.init + 90000;
defaults.profile = 'compact';
defaults.maxCandidates = 2;
defaults.randomizeW = false;

names = fieldnames(defaults);
for ii = 1:numel(names)
    name = names{ii};
    if ~isfield(opts,name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
opts.maxCandidates = max(1,round(opts.maxCandidates));
end
