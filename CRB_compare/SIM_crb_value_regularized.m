function [JcrbReg, info] = SIM_crb_value_regularized(ch, P, Rx, alg)
%SIM_CRB_VALUE_REGULARIZED Evaluate a numerically stabilized CRB surrogate.
%
% This function is for internal search/diagnostic plotting only. It does not
% replace SIM_crb_value, which remains the exact unregularized physical CRB
% used for retained data and final acceptance checks.
%
% Areg = A + deltaA*I,  deltaA = regRel*meanEig(A)
% Breg = B + deltaB*I,  deltaB = regRel*meanEig(B)
%
% where A = P*G*Rx*G^H*P^H and B = P*G*G^H*P^H.

if nargin < 4 || isempty(alg)
    alg = struct();
end
regRel = get_option(alg,'crbInternalRegRel',1e-9);

JcrbReg = Inf;
info = initialize_info();
info.regRel = regRel;

Cx = ch.G * Rx * ch.G';
C0 = ch.G * ch.G';
Cx = (Cx + Cx')/2;
C0 = (C0 + C0')/2;
A = P * Cx * P';
B = P * C0 * P';
A = (A + A')/2;
B = (B + B')/2;

info.Cx = Cx;
info.C0 = C0;
info.A = A;
info.B = B;

if any(~isfinite(A(:))) || any(~isfinite(B(:)))
    info.invalidReason = 'A or B contains a nonfinite entry.';
    return;
end

[UA,dA] = eig(A,'vector');
[UB,dB] = eig(B,'vector');
dA = real(dA(:));
dB = real(dB(:));
[dA,ordA] = sort(dA,'descend');
[dB,ordB] = sort(dB,'descend');
UA = UA(:,ordA);
UB = UB(:,ordB);

maxEigA = max(dA);
maxEigB = max(dB);
if ~(isfinite(maxEigA) && maxEigA > 0 && ...
        isfinite(maxEigB) && maxEigB > 0)
    info.invalidReason = 'A or B has no positive largest eigenvalue.';
    return;
end

nA = numel(dA);
nB = numel(dB);
scaleA = max(real(trace(A))/max(nA,1),maxEigA*eps);
scaleB = max(real(trace(B))/max(nB,1),maxEigB*eps);

if ~(isfinite(regRel) && regRel >= 0)
    info.invalidReason = 'crbInternalRegRel must be finite and nonnegative.';
    return;
end

deltaA = regRel*scaleA;
deltaB = regRel*scaleB;
dAreg = dA + deltaA;
dBreg = dB + deltaB;

info.regAbsA = deltaA;
info.regAbsB = deltaB;
info.minEigRawA = min(dA);
info.minEigRawB = min(dB);
info.maxEigRawA = maxEigA;
info.maxEigRawB = maxEigB;
info.minEigRegA = min(dAreg);
info.minEigRegB = min(dBreg);

if any(~isfinite(dAreg)) || any(~isfinite(dBreg)) || ...
        any(dAreg <= 0) || any(dBreg <= 0)
    info.invalidReason = ['Regularized spectrum is not positive. ', ...
        'The raw information matrix is too indefinite for the configured loading.'];
    return;
end

[Ainv,traceAinv,logTraceAinv,okA] = ...
    inverse_from_spectrum(UA,dAreg);
[Binv,traceBinv,logTraceBinv,okB] = ...
    inverse_from_spectrum(UB,dBreg);
if ~(okA && okB)
    info.invalidReason = 'Regularized inverse trace is nonfinite.';
    return;
end

logJ = logTraceAinv + logTraceBinv;
if ~(isfinite(logJ) && logJ < log(realmax))
    info.invalidReason = 'Regularized CRB logarithm is nonfinite or overflowed.';
    return;
end

Jtmp = exp(logJ);
if ~(isfinite(Jtmp) && Jtmp > 0)
    info.invalidReason = 'Regularized CRB is nonfinite or nonpositive.';
    return;
end

JcrbReg = Jtmp;
info.Ainv = Ainv;
info.Binv = Binv;
info.traceAinv = traceAinv;
info.traceBinv = traceBinv;
info.logTraceAinv = logTraceAinv;
info.logTraceBinv = logTraceBinv;
info.logCRB = logJ;
info.valid = true;
info.invalidReason = '';
end

function [Xinv,traceInv,logTraceInv,ok] = inverse_from_spectrum(U,d)
d = real(d(:));
ok = all(isfinite(d)) && all(d > 0);
Xinv = [];
traceInv = Inf;
logTraceInv = Inf;
if ~ok
    return;
end

dmax = max(d);
mu = d/max(dmax,realmin);
invMu = 1./mu;
traceInv = sum(invMu)/dmax;
logTraceInv = log(sum(invMu)) - log(dmax);
invEig = invMu/dmax;
Xinv = U*diag(invEig)*U';
Xinv = (Xinv+Xinv')/2;
ok = all(isfinite(Xinv(:))) && ...
    isfinite(traceInv) && traceInv > 0 && isfinite(logTraceInv);
end

function info = initialize_info()
info = struct();
info.valid = false;
info.invalidReason = 'Not evaluated.';
info.regRel = NaN;
info.regAbsA = NaN;
info.regAbsB = NaN;
info.Cx = [];
info.C0 = [];
info.A = [];
info.B = [];
info.Ainv = [];
info.Binv = [];
info.traceAinv = Inf;
info.traceBinv = Inf;
info.logTraceAinv = Inf;
info.logTraceBinv = Inf;
info.logCRB = Inf;
info.minEigRawA = NaN;
info.minEigRawB = NaN;
info.maxEigRawA = NaN;
info.maxEigRawB = NaN;
info.minEigRegA = NaN;
info.minEigRegB = NaN;
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
