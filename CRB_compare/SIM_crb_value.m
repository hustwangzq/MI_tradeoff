function [Jcrb, info] = SIM_crb_value(ch, P, Rx, alg)
%SIM_CRB_VALUE Evaluate the SIM CRB objective without regularization.
%
%   J = tr(A^{-1}) tr(B^{-1}),
%   A = P G Rx G^H P^H,
%   B = P G G^H P^H.
%
% The inverse traces are evaluated through the Hermitian eigendecomposition.
% Each eigenvalue set is scaled by its largest eigenvalue before reciprocals
% are formed. This is algebraically identical to the original CRB and does
% not clip, load, or otherwise modify any eigenvalue.
%
% Two relative-eigenvalue thresholds are supported:
%   crbWarnMinRelEig : diagnostic warning only;
%   crbHardMinRelEig : below this value the CRB is treated as numerically
%                      unreliable and the candidate is rejected.

if nargin < 4 || isempty(alg)
    alg = struct();
end
hardTol = get_option(alg,'crbHardMinRelEig', ...
    get_option(alg,'crbMinRcond',1e-12));
warnTol = get_option(alg,'crbWarnMinRelEig',max(100*hardTol,1e-10));

Jcrb = Inf;
info = initialize_info();
info.hardThreshold = hardTol;
info.warningThreshold = warnTol;

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

info.minEigA = min(dA);
info.maxEigA = max(dA);
info.minEigB = min(dB);
info.maxEigB = max(dB);
info.rcondA = rcond(A);
info.rcondB = rcond(B);

if ~(isfinite(info.maxEigA) && info.maxEigA > 0 && ...
        isfinite(info.maxEigB) && info.maxEigB > 0)
    info.invalidReason = 'A or B has no positive largest eigenvalue.';
    return;
end

muA = dA/info.maxEigA;
muB = dB/info.maxEigB;
info.normalizedEigA = muA;
info.normalizedEigB = muB;
info.relMinEigA = min(muA);
info.relMinEigB = min(muB);
info.conditionWarning = info.relMinEigA < warnTol || ...
    info.relMinEigB < warnTol;

if any(~isfinite(muA)) || any(~isfinite(muB)) || ...
        info.minEigA <= 0 || info.minEigB <= 0
    info.invalidReason = 'A or B is not positive definite.';
    return;
end
if info.relMinEigA < hardTol || info.relMinEigB < hardTol
    info.invalidReason = sprintf([ ...
        'Relative minimum eigenvalue below hard threshold: ', ...
        'A %.3e, B %.3e, threshold %.3e.'], ...
        info.relMinEigA,info.relMinEigB,hardTol);
    return;
end

% Exact spectral inverses. The reciprocal is formed in normalized units and
% then restored by the largest-eigenvalue scale.
invEigA = (1./muA)/info.maxEigA;
invEigB = (1./muB)/info.maxEigB;
Ainv = UA*diag(invEigA)*UA';
Binv = UB*diag(invEigB)*UB';
Ainv = (Ainv + Ainv')/2;
Binv = (Binv + Binv')/2;

traceAinv = sum(invEigA);
traceBinv = sum(invEigB);
logTraceAinv = log(sum(1./muA)) - log(info.maxEigA);
logTraceBinv = log(sum(1./muB)) - log(info.maxEigB);
logJ = logTraceAinv + logTraceBinv;

info.Ainv = Ainv;
info.Binv = Binv;
% Scale-free inverse shapes used by the phase gradient. These avoid forming
% products of two very large physical inverses for strongly attenuated
% multi-hop rRIS channels. The resulting log-gradient is algebraically
% identical to the physical expression.
info.AinvNormalized = info.maxEigA*Ainv;
info.BinvNormalized = info.maxEigB*Binv;
info.traceAinvNormalized = sum(1./muA);
info.traceBinvNormalized = sum(1./muB);
info.traceAinv = traceAinv;
info.traceBinv = traceBinv;
info.logTraceAinv = logTraceAinv;
info.logTraceBinv = logTraceBinv;
info.logCRB = logJ;

if any(~isfinite(Ainv(:))) || any(~isfinite(Binv(:))) || ...
        ~(isfinite(traceAinv) && traceAinv > 0) || ...
        ~(isfinite(traceBinv) && traceBinv > 0) || ...
        ~(isfinite(logJ) && logJ < log(realmax))
    info.invalidReason = 'The exact inverse trace is nonfinite or overflowed.';
    return;
end

Jtmp = exp(logJ);
if ~(isfinite(Jtmp) && Jtmp > 0)
    info.invalidReason = 'The CRB product is nonfinite or nonpositive.';
    return;
end

Jcrb = Jtmp;
info.valid = true;
info.invalidReason = '';
end

function info = initialize_info()
info = struct();
info.valid = false;
info.invalidReason = 'Not evaluated.';
info.conditionWarning = false;
info.hardThreshold = NaN;
info.warningThreshold = NaN;
info.Cx = [];
info.C0 = [];
info.A = [];
info.B = [];
info.Ainv = [];
info.Binv = [];
info.AinvNormalized = [];
info.BinvNormalized = [];
info.traceAinvNormalized = Inf;
info.traceBinvNormalized = Inf;
info.traceAinv = Inf;
info.traceBinv = Inf;
info.logTraceAinv = Inf;
info.logTraceBinv = Inf;
info.logCRB = Inf;
info.rcondA = NaN;
info.rcondB = NaN;
info.minEigA = NaN;
info.maxEigA = NaN;
info.relMinEigA = NaN;
info.minEigB = NaN;
info.maxEigB = NaN;
info.relMinEigB = NaN;
info.normalizedEigA = [];
info.normalizedEigB = [];
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
