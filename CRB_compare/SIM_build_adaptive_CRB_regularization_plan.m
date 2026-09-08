function plan = SIM_build_adaptive_CRB_regularization_plan(ch,P,alg)
%SIM_BUILD_ADAPTIVE_CRB_REGULARIZATION_PLAN Build numerical fallback plan.
%
% The plan changes only internal search regularization and the optional
% sensing excitation floor. Every returned covariance is still certified
% with the exact unregularized physical CRB by SIM_solve_CRB_SINR_SDR.
%
% Adaptive policy:
%   1) start from no A-side loading and a B loading inferred from raw B
%      conditioning;
%   2) gradually increase A/B numerical loading;
%   3) simultaneously relax the artificial sensing-eigenvalue floor so the
%      search can approach the pure-communication SINR boundary;
%   4) use a fixed regularization fallback only as the final attempt.

if nargin < 3 || isempty(alg)
    alg = struct('robustSDRFixedFallbackRegRel',1e-7);
end

N = size(P,1);
E = P*ch.G;
B = E*E';
B = (B+B')/2;
eigB = real(eig(B));
maxEigB = max(eigB);
minEigB = min(eigB);
if isfinite(maxEigB) && maxEigB > 0
    relB = minEigB/maxEigB;
else
    relB = NaN;
end

fixedReg = get_option(alg,'robustSDRFixedFallbackRegRel',1e-7);
fixedReg = sanitize_positive(fixedReg,1e-7);
baseFloor = get_option(alg,'sensingEigFloorRel',0);
baseFloor = sanitize_nonnegative(baseFloor,0);

% B-side adaptive starting value. This is intentionally mild: B must still
% pass the RAW physical-rank check in the base solver.
if ~isfinite(relB)
    adaptiveB0 = 1e-8;
elseif relB >= 1e-5
    adaptiveB0 = 0;
elseif relB >= 1e-7
    adaptiveB0 = 1e-10;
elseif relB >= 1e-9
    adaptiveB0 = 1e-9;
else
    adaptiveB0 = 1e-8;
end

adaptiveA = [0,1e-10,1e-9,1e-8,fixedReg];
adaptiveB = [adaptiveB0,max(adaptiveB0,1e-10), ...
    max(adaptiveB0,1e-9),max(adaptiveB0,1e-8),fixedReg];

% The floor is NOT part of the reported CRB definition. It is progressively
% relaxed so that high-SINR points are not artificially capped by a fixed
% sensing-covariance floor.
floorList = [baseFloor,min(baseFloor,1e-10),0,0,0];

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

basePrecision=lower(char(string(get_option( ...
    alg,'robustSDRPrecisionMode','default'))));
if ~ismember(basePrecision,{'default','high','best'}), basePrecision='default'; end
precisionList = repmat({basePrecision},1,numel(adaptiveA));
if logical(get_option(alg,'robustSDRFinalBestRetry',false))
    adaptiveA(end+1) = fixedReg;
    adaptiveB(end+1) = fixedReg;
    floorList(end+1) = 0;
    finalPrecision=lower(char(string(get_option( ...
        alg,'robustSDRFinalPrecisionMode','best'))));
    if ~ismember(finalPrecision,{'default','high','best'}), finalPrecision='best'; end
    precisionList{end+1} = finalPrecision;
end

plan = struct( ...
    'numAttempts',numel(adaptiveA), ...
    'regRelA',adaptiveA, ...
    'regRelB',adaptiveB, ...
    'sensingFloorRel',floorList, ...
    'precision',{precisionList}, ...
    'schurScale',scaleValue*ones(1,numel(adaptiveA)), ...
    'rawRelMinEigB',relB, ...
    'fixedFallbackRegRel',fixedReg, ...
    'mode','adaptive-A/B-with-floor-relaxation');
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end

function value = sanitize_positive(value,defaultValue)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value > 0)
    value = defaultValue;
end
end

function value = sanitize_nonnegative(value,defaultValue)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    value = defaultValue;
end
end
