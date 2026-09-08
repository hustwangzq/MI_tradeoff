function sdr = SIM_solve_CRB_SINR_SDR(params, ch, P, gammaTarget, alg)
%SIM_SOLVE_CRB_SINR_SDR Fixed-phase CRB-SINR covariance SDR.
%
% Reported/accepted CRB is ALWAYS the exact unregularized physical CRB.
% Internal regularization is used only to stabilize the search model:
%   - B-side loading stabilizes sensing-channel whitening;
%   - A-side loading stabilizes the inverse Schur block in whitened space.
% Neither loading modifies the communication constraints nor the final raw
% CRB post-check.
%
% The Schur block supports numerical scaling through alg.cvxSchurScale.

if nargin < 5 || isempty(alg)
    alg = struct('cvxPrecisionMode','high');
end

sdr = initialize_output();

%% Basic checks and defaults
requiredParams = {'K','Nt','P0'};
for ii = 1:numel(requiredParams)
    if ~isfield(params,requiredParams{ii})
        sdr.status = sprintf('Missing params.%s.',requiredParams{ii});
        return;
    end
end
if ~isfield(ch,'G') || ~isfield(ch,'hUsers') || ~isfield(ch,'sigma_c2')
    sdr.status = 'Channel structure must contain G, hUsers, and sigma_c2.';
    return;
end

K = params.K;
Nt = params.Nt;
N = size(P,1);
P0 = params.P0;

gammaTarget = gammaTarget(:);
if numel(gammaTarget) ~= K || any(~isfinite(gammaTarget)) || any(gammaTarget < 0)
    sdr.status = 'gammaTarget must contain K finite nonnegative entries.';
    return;
end
if size(ch.G,2) ~= Nt
    sdr.status = sprintf('size(ch.G,2)=%d is inconsistent with Nt=%d.', ...
        size(ch.G,2),Nt);
    return;
end
if size(P,2) ~= size(ch.G,1)
    sdr.status = 'P and ch.G have incompatible dimensions.';
    return;
end
if ~(isfinite(P0) && P0 > 0)
    sdr.status = 'params.P0 must be finite and positive.';
    return;
end

crbHardMinRelEig = get_option(alg,'crbHardMinRelEig', ...
    get_option(alg,'crbMinRcond',1e-12));
crbWarnMinRelEig = get_option(alg,'crbWarnMinRelEig', ...
    max(100*crbHardMinRelEig,1e-10));
sensingEigFloorRel = get_option(alg,'sensingEigFloorRel',0);

% Backward compatible common switch plus separate A/B internal loadings.
useInternalCRBRegularization = logical(get_option( ...
    alg,'useInternalCRBRegularization',false));
legacyRegRel = get_option(alg,'crbInternalRegRel',1e-9);
crbInternalRegRelA = get_option(alg,'crbInternalRegRelA',legacyRegRel);
crbInternalRegRelB = get_option(alg,'crbInternalRegRelB',legacyRegRel);
if ~useInternalCRBRegularization
    crbInternalRegRelA = 0;
    crbInternalRegRelB = 0;
end
crbInternalRegRelA = sanitize_nonnegative_scalar(crbInternalRegRelA,0);
crbInternalRegRelB = sanitize_nonnegative_scalar(crbInternalRegRelB,0);

cvxQuiet = logical(get_option(alg,'cvxQuiet',true));
cvxPowerTol = get_option(alg,'cvxPowerTol',1e-8);
postSINRTol = get_option(alg,'sdrPostSINRTol',1e-5);
useMixedSINRTolerance = ...
    (isfield(alg,'sdrPostSINRRelTol') && ~isempty(alg.sdrPostSINRRelTol)) || ...
    (isfield(alg,'sdrPostSINRAbsTol') && ~isempty(alg.sdrPostSINRAbsTol));
postSINRRelTol = get_option(alg,'sdrPostSINRRelTol',1e-5);
postSINRAbsTol = get_option(alg,'sdrPostSINRAbsTol',1e-10);
postPSDTol = get_option(alg,'sdrPostPSDTol',1e-7);

validateattributes(postSINRTol,{'numeric'}, ...
    {'real','finite','scalar','nonnegative'},mfilename,'sdrPostSINRTol');
validateattributes(postSINRRelTol,{'numeric'}, ...
    {'real','finite','scalar','nonnegative'},mfilename,'sdrPostSINRRelTol');
validateattributes(postSINRAbsTol,{'numeric'}, ...
    {'real','finite','scalar','nonnegative'},mfilename,'sdrPostSINRAbsTol');

if isfield(alg,'cvxPrecisionMode') && ~isempty(alg.cvxPrecisionMode)
    cvxPrecisionMode = lower(char(string(alg.cvxPrecisionMode)));
elseif logical(get_option(alg,'cvxHighPrecision',true))
    cvxPrecisionMode = 'high';
else
    cvxPrecisionMode = 'default';
end
if ~ismember(cvxPrecisionMode,{'default','high','best'})
    error('SIM_solve_CRB_SINR_SDR:BadPrecision', ...
        'alg.cvxPrecisionMode must be default, high, or best.');
end
cvxSchurScale = get_option(alg,'cvxSchurScale',1);
validateattributes(cvxSchurScale,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'alg.cvxSchurScale');

%% Fixed sensing channel and B-side whitening
E = P * ch.G;
B = E * E';
B = (B + B')/2;
if any(~isfinite(B(:)))
    sdr.status = 'Nonfinite fixed sensing Gram matrix B.';
    return;
end

[U,Dmat] = eig(B,'vector');
d = real(Dmat(:));
[d,ord] = sort(d,'descend');
U = U(:,ord);
maxEigB = max(d);
minEigB = min(d);
if ~(isfinite(maxEigB) && maxEigB > 0)
    sdr.status = 'Fixed sensing Gram matrix B has no positive eigenvalue.';
    return;
end
relativeMinEigB = minEigB/maxEigB;
sdr.maxEigB = maxEigB;
sdr.minEigB = minEigB;
sdr.relMinEigB = relativeMinEigB;
sdr.rcondB = rcond(B);
sdr.rankB = sum(d > crbHardMinRelEig*maxEigB);
sdr.conditionWarningB = relativeMinEigB < crbWarnMinRelEig;

% The raw B must remain physically invertible. Internal loading cannot turn
% a truly invalid physical point into an accepted result.
if minEigB <= 0 || ~isfinite(relativeMinEigB) || ...
        relativeMinEigB < crbHardMinRelEig
    sdr.status = sprintf(['Fixed sensing Gram matrix is singular or below ', ...
        'the physical relative-eigenvalue threshold: relMinEig(B)=%.3e, ', ...
        'threshold=%.3e.'],relativeMinEigB,crbHardMinRelEig);
    return;
end

scaleB = max(real(trace(B))/max(N,1),maxEigB*eps);
regAbsB = crbInternalRegRelB*scaleB;
dWhite = d + regAbsB;
maxEigBWhite = max(dWhite);
if any(~isfinite(dWhite)) || any(dWhite <= 0) || ...
        ~(isfinite(maxEigBWhite) && maxEigBWhite > 0)
    sdr.status = 'Invalid internally regularized sensing spectrum.';
    return;
end

muBWhite = dWhite/maxEigBWhite;
invSqrtEigB = (1./sqrt(muBWhite))/sqrt(maxEigBWhite);
invEigB = (1./muBWhite)/maxEigBWhite;
BinvHalf = diag(invSqrtEigB) * U';
F = full(BinvHalf * E);

Binv = U * diag(invEigB) * U';
Binv = (Binv + Binv')/2;
traceBinv = sum(invEigB);
if ~(isfinite(traceBinv) && traceBinv > 0)
    sdr.status = 'Invalid internal trace(B^{-1}) during sensing whitening.';
    return;
end
WB = Binv/traceBinv;
WB = (WB + WB')/2;

whiteningError = norm(F*F' - eye(N),'fro')/max(N,1);
sdr.whiteningError = whiteningError;
sdr.traceBinv = traceBinv;
sdr.Fwhite = F;
sdr.inverseWeight = WB;
sdr.internalRegularizationUsed = (crbInternalRegRelA > 0) || (regAbsB > 0);
sdr.internalRegRel = max(crbInternalRegRelA,crbInternalRegRelB);
sdr.internalRegRelA = crbInternalRegRelA;
sdr.internalRegRelB = crbInternalRegRelB;
sdr.internalRegAbsAWhite = crbInternalRegRelA;
sdr.internalRegAbsB = regAbsB;
sdr.sensingEigFloorRelUsed = sensingEigFloorRel;

%% Per-user communication scaling
Hscaled = zeros(Nt,Nt,K);
noiseScaled = zeros(K,1);
gUsers = zeros(Nt,K);
for k = 1:K
    if k > numel(ch.hUsers) || isempty(ch.hUsers{k})
        sdr.status = sprintf('Missing ch.hUsers{%d}.',k);
        return;
    end
    gk = ch.G' * P' * ch.hUsers{k};
    gk = gk(:);
    if numel(gk) ~= Nt || any(~isfinite(gk))
        sdr.status = sprintf('Invalid effective communication channel for user %d.',k);
        return;
    end
    gUsers(:,k) = gk;
    Hk = gk*gk';
    Hk = (Hk + Hk')/2;
    userScale = P0*real(trace(Hk)) + ch.sigma_c2;
    userScale = max(userScale,realmin);
    Hscaled(:,:,k) = (P0/userScale)*Hk;
    noiseScaled(k) = ch.sigma_c2/userScale;
end

%% CVX problem
try
    if cvxQuiet
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end

        if isfield(alg,'cvxSolver') && ~isempty(alg.cvxSolver)
            cvx_solver(char(string(alg.cvxSolver)));
        end
        switch cvxPrecisionMode
            case 'default'
                cvx_precision default
            case 'high'
                cvx_precision high
            case 'best'
                cvx_precision best
        end

        variable Xcomm(Nt,Nt,K) hermitian semidefinite
        variable Xs(Nt,Nt) hermitian semidefinite
        variable Zscaled(N,N) hermitian semidefinite

        expression Xtotal(Nt,Nt)
        expression Cwhite(N,N)
        expression CwhiteReg(N,N)

        Xtotal = Xs;
        for k = 1:K
            Xtotal = Xtotal + Xcomm(:,:,k);
        end

        Cwhite = F*Xtotal*F';
        Cwhite = 0.5*(Cwhite + Cwhite');

        % A-side adaptive regularization is applied ONLY to the inverse
        % search model. In whitened/normalized coordinates Cwhite has O(1)
        % natural scale, so crbInternalRegRelA is dimensionless here.
        CwhiteReg = Cwhite + crbInternalRegRelA*eye(N);

        minimize( cvxSchurScale*real(trace(WB*Zscaled)) )

        subject to
            real(trace(Xtotal)) <= 1;

            [cvxSchurScale*CwhiteReg, eye(N); eye(N), Zscaled] >= 0;

            % Optional excitation floor is a feasibility restriction. Robust
            % continuation is allowed to relax it adaptively, but the final
            % physical CRB is still checked without any loading.
            if sensingEigFloorRel > 0
                Cwhite >= sensingEigFloorRel*eye(N);
            end

            for k = 1:K
                real(trace(Hscaled(:,:,k)*Xcomm(:,:,k))) >= ...
                    gammaTarget(k) * ( ...
                    real(trace(Hscaled(:,:,k)*(Xtotal-Xcomm(:,:,k)))) ...
                    + noiseScaled(k) );
            end
    cvx_end

catch ME
    try
        cvx_clear;
    catch
    end
    sdr.status = sprintf('CVX exception: %s',ME.message);
    sdr.cvxStatus = 'Exception';
    return;
end

sdr.cvxStatus = char(string(cvx_status));
sdr.cvxOptval = cvx_optval;
sdr.cvxPrecisionMode = cvxPrecisionMode;
sdr.cvxSchurScale = cvxSchurScale;
if exist('cvx_slvtol','var') && isscalar(cvx_slvtol) && isfinite(cvx_slvtol)
    sdr.cvxSlvtol = cvx_slvtol;
else
    sdr.cvxSlvtol = NaN;
end

cvxSolved = contains(lower(sdr.cvxStatus),'solved');
if ~cvxSolved
    sdr.status = sdr.cvxStatus;
    return;
end

%% Recover covariances in the original power scale
Qcomm = zeros(Nt,Nt,K);
for k = 1:K
    Qcomm(:,:,k) = P0*full(Xcomm(:,:,k));
    Qcomm(:,:,k) = (Qcomm(:,:,k)+Qcomm(:,:,k)')/2;
end
Rs = P0*full(Xs);
Rs = (Rs+Rs')/2;
Rx = Rs;
for k = 1:K
    Rx = Rx + Qcomm(:,:,k);
end
Rx = (Rx+Rx')/2;

%% Exact numerical post-check in ORIGINAL variables
powerUsed = real(trace(Rx));
sinr = zeros(K,1);
sinrResidual = zeros(K,1);
for k = 1:K
    gk = gUsers(:,k);
    desired = real(gk'*Qcomm(:,:,k)*gk);
    interference = real(gk'*Rs*gk);
    for j = 1:K
        if j ~= k
            interference = interference + real(gk'*Qcomm(:,:,j)*gk);
        end
    end
    denom = interference + ch.sigma_c2;
    sinr(k) = desired/max(denom,realmin);
    sinrResidual(k) = desired - gammaTarget(k)*denom;
end

minEigQ = Inf(K,1);
for k = 1:K
    minEigQ(k) = min(real(eig(Qcomm(:,:,k))));
end
minEigRs = min(real(eig(Rs)));

% IMPORTANT: raw/unregularized CRB only.
[Jcrb,crbInfo] = SIM_crb_value(ch,P,Rx,alg);

powerTolAbs = max(cvxPowerTol,1e-7*P0);
if useMixedSINRTolerance
    sinrTolerance = postSINRAbsTol + postSINRRelTol.*abs(gammaTarget);
    sinrToleranceMode = 'mixed';
else
    sinrScale = max(abs(gammaTarget),1);
    sinrTolerance = postSINRTol.*sinrScale;
    sinrToleranceMode = 'legacy';
end
sinrFeasible = all(sinr + sinrTolerance >= gammaTarget);
powerFeasible = powerUsed <= P0 + powerTolAbs;
psdScale = max([1,norm(Rx,2),P0]);
psdFeasible = all(minEigQ >= -postPSDTol*psdScale) && ...
    minEigRs >= -postPSDTol*psdScale;
crbFeasible = crbInfo.valid && isfinite(Jcrb);

sdr.Qcomm = Qcomm;
sdr.Rs = Rs;
sdr.Rx = Rx;
sdr.sinr = sinr;
sdr.sinrResidual = sinrResidual;
sdr.sinrTolerance = sinrTolerance;
sdr.sinrToleranceMode = sinrToleranceMode;
sdr.postSINRTol = postSINRTol;
sdr.postSINRRelTol = postSINRRelTol;
sdr.postSINRAbsTol = postSINRAbsTol;
sdr.power = powerUsed;
sdr.minEigQcomm = minEigQ;
sdr.minEigRs = minEigRs;
sdr.CRB = Jcrb;
sdr.crbInfo = crbInfo;
sdr.objectiveWhitened = real(cvx_optval);
sdr.objectiveCRBEstimate = real(cvx_optval)*(traceBinv^2/P0);
sdr.powerFeasible = powerFeasible;
sdr.sinrFeasible = sinrFeasible;
sdr.psdFeasible = psdFeasible;
sdr.crbFeasible = crbFeasible;
sdr.success = cvxSolved && powerFeasible && sinrFeasible ...
    && psdFeasible && crbFeasible;

if sdr.success
    sdr.status = sdr.cvxStatus;
else
    reasons = cell(1,0);
    if ~powerFeasible
        reasons{end+1} = sprintf('power %.3e > %.3e',powerUsed,P0); %#ok<AGROW>
    end
    if ~sinrFeasible
        reasons{end+1} = sprintf('min achieved SINR %.3e below target %.3e', ...
            min(sinr),min(gammaTarget)); %#ok<AGROW>
    end
    if ~psdFeasible
        reasons{end+1} = sprintf('PSD check failed: minEigQ=%.3e, minEigRs=%.3e', ...
            min(minEigQ),minEigRs); %#ok<AGROW>
    end
    if ~crbFeasible
        reasons{end+1} = sprintf('raw CRB invalid: rcondA=%.3e, rcondB=%.3e', ...
            safe_info(crbInfo,'rcondA'),safe_info(crbInfo,'rcondB')); %#ok<AGROW>
    end
    if isempty(reasons)
        sdr.status = ['CVX reported ',sdr.cvxStatus, ...
            ', but the exact post-check did not accept the solution.'];
    else
        sdr.status = ['Post-check failed after ',sdr.cvxStatus,': ', ...
            strjoin(reasons,'; ')];
    end
end
end

function sdr = initialize_output()
% Fully populated fixed-schema scalar structure. No empty struct assignment.
crbInfo0 = struct('valid',false,'rcondA',NaN,'rcondB',NaN, ...
    'relMinEigA',NaN,'relMinEigB',NaN);
sdr = struct( ...
    'success',false, ...
    'status','Not solved', ...
    'cvxStatus','Not run', ...
    'cvxOptval',NaN, ...
    'cvxSlvtol',NaN, ...
    'cvxPrecisionMode','not-set', ...
    'cvxSchurScale',NaN, ...
    'Qcomm',[], ...
    'Rs',[], ...
    'Rx',[], ...
    'sinr',[], ...
    'sinrResidual',[], ...
    'sinrTolerance',[], ...
    'sinrToleranceMode','not-set', ...
    'postSINRTol',NaN, ...
    'postSINRRelTol',NaN, ...
    'postSINRAbsTol',NaN, ...
    'power',NaN, ...
    'CRB',Inf, ...
    'crbInfo',crbInfo0, ...
    'objectiveWhitened',NaN, ...
    'objectiveCRBEstimate',NaN, ...
    'powerFeasible',false, ...
    'sinrFeasible',false, ...
    'psdFeasible',false, ...
    'crbFeasible',false, ...
    'minEigQcomm',[], ...
    'minEigRs',NaN, ...
    'maxEigB',NaN, ...
    'minEigB',NaN, ...
    'relMinEigB',NaN, ...
    'rcondB',NaN, ...
    'rankB',NaN, ...
    'whiteningError',NaN, ...
    'traceBinv',NaN, ...
    'Fwhite',[], ...
    'inverseWeight',[], ...
    'conditionWarningB',false, ...
    'internalRegularizationUsed',false, ...
    'internalRegRel',0, ...
    'internalRegRelA',0, ...
    'internalRegRelB',0, ...
    'internalRegAbsAWhite',0, ...
    'internalRegAbsB',0, ...
    'sensingEigFloorRelUsed',NaN);
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end

function value = sanitize_nonnegative_scalar(value,defaultValue)
if ~(isnumeric(value) && isscalar(value) && isfinite(value) && value >= 0)
    value = defaultValue;
end
end

function value = safe_info(info,name)
value = NaN;
if isstruct(info) && isfield(info,name) && isscalar(info.(name))
    value = info.(name);
end
end
