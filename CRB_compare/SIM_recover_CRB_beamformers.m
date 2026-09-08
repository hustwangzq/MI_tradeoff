function rec = SIM_recover_CRB_beamformers(params, ch, P, Qcomm, RsSDR, alg)
%SIM_RECOVER_CRB_BEAMFORMERS Covariance-preserving rank-one recovery.
%
% For user k,
%   w_k = Q_k g_k / sqrt(g_k^H Q_k g_k).
% The residual Q_k-w_k w_k^H is transferred to the sensing covariance.
% Under exact arithmetic this preserves the total covariance and all user
% SINRs. This implementation deliberately does NOT fall back to a principal
% eigenvector when the denominator becomes numerically unreliable; such an
% SDR point is reported as an invalid recovery so the AO layer can reject it.

Nt = params.Nt;
K = params.K;

softTol = get_option(alg,'recoveryPSDSoftTol',1e-10);
hardTol = get_option(alg,'recoveryPSDHardTol',1e-7);
denomTol = get_option(alg,'recoveryDenomRelTol', ...
    get_option(alg,'recoveryEigTol',1e-12));
covGapTol = get_option(alg,'recoveryCovarianceGapTol',1e-7);
desiredGapTol = get_option(alg,'recoveryDesiredGapTol',1e-8);
leakTol = get_option(alg,'recoverySelfLeakTol',1e-8);

rec = make_empty_recovery(Nt,K);

if any(~isfinite(Qcomm(:))) || any(~isfinite(RsSDR(:)))
    rec.status = 'nonfinite-SDR-covariance';
    return;
end

RxSDR = (RsSDR+RsSDR')/2;
for k = 1:K
    RxSDR = RxSDR + (Qcomm(:,:,k)+Qcomm(:,:,k)')/2;
end
RxSDR = (RxSDR+RxSDR')/2;

Wc = zeros(Nt,K);
RsRecovered = (RsSDR+RsSDR')/2;
residuals = zeros(Nt,Nt,K);
minResidualRelEig = Inf;
maxResidualCorrectionRel = 0;
maxDesiredSignalGapRel = 0;
maxSelfResidualLeakRel = 0;

for k = 1:K
    Qk = (Qcomm(:,:,k)+Qcomm(:,:,k)')/2;
    gk = ch.G' * P' * ch.hUsers{k};
    if any(~isfinite(Qk(:))) || any(~isfinite(gk(:)))
        rec.status = sprintf('nonfinite-user-data-k%d',k);
        return;
    end

    denom = real(gk'*Qk*gk);
    denomScale = max(real(trace(Qk))*norm(gk)^2,realmin);
    if ~isfinite(denom) || denom <= denomTol*denomScale
        rec.status = sprintf('recovery-denominator-too-small-k%d',k);
        rec.failedUser = k;
        rec.recoveryDenominator = denom;
        rec.recoveryDenominatorScale = denomScale;
        return;
    end

    wk = Qk*gk/sqrt(denom);
    if any(~isfinite(wk))
        rec.status = sprintf('nonfinite-recovered-beam-k%d',k);
        rec.failedUser = k;
        return;
    end
    Wc(:,k) = wk;

    desiredRecovered = abs(gk'*wk)^2;
    desiredGapRel = abs(desiredRecovered-denom)/max(abs(denom),realmin);
    maxDesiredSignalGapRel = max(maxDesiredSignalGapRel,desiredGapRel);

    residualRaw = (Qk-wk*wk' + (Qk-wk*wk')')/2;
    [residual,psdInfo] = sanitize_psd(residualRaw,softTol,hardTol);
    minResidualRelEig = min(minResidualRelEig,psdInfo.minRelEigBefore);
    maxResidualCorrectionRel = max(maxResidualCorrectionRel, ...
        psdInfo.correctionRel);
    if ~psdInfo.valid
        rec.status = sprintf('materially-indefinite-residual-k%d',k);
        rec.failedUser = k;
        rec.minResidualRelEigBefore = psdInfo.minRelEigBefore;
        return;
    end

    selfLeak = abs(real(gk'*residual*gk))/max(abs(denom),realmin);
    maxSelfResidualLeakRel = max(maxSelfResidualLeakRel,selfLeak);
    residuals(:,:,k) = residual;
    RsRecovered = RsRecovered + residual;
end

[RsRecovered,rsInfo] = sanitize_psd( ...
    (RsRecovered+RsRecovered')/2,softTol,hardTol);
if ~rsInfo.valid
    rec.status = 'materially-indefinite-recovered-sensing-covariance';
    rec.minRsRelEigBefore = rsInfo.minRelEigBefore;
    return;
end

[U,D] = eig((RsRecovered+RsRecovered')/2);
d = real(diag(D));
if any(~isfinite(d))
    rec.status = 'nonfinite-sensing-eigendecomposition';
    return;
end
d = max(d,0);
Wr = U*diag(sqrt(d));

W = [Wc,Wr];
RxRecovered = (W*W' + (W*W')')/2;

covarianceGapRel = norm(RxRecovered-RxSDR,'fro') / ...
    max(norm(RxSDR,'fro'),realmin);
powerSDR = real(trace(RxSDR));
powerRecovered = real(trace(RxRecovered));
powerGapRel = abs(powerRecovered-powerSDR)/max(abs(powerSDR),realmin);

rec.W = W;
rec.Wc = Wc;
rec.Wr = Wr;
rec.RsRecovered = RsRecovered;
rec.RxRecovered = RxRecovered;
rec.RxSDR = RxSDR;
rec.residuals = residuals;
rec.power = powerRecovered;
rec.powerSDR = powerSDR;
rec.covarianceGapRel = covarianceGapRel;
rec.powerGapRel = powerGapRel;
rec.maxDesiredSignalGapRel = maxDesiredSignalGapRel;
rec.maxSelfResidualLeakRel = maxSelfResidualLeakRel;
rec.minResidualRelEigBefore = minResidualRelEig;
rec.maxResidualCorrectionRel = maxResidualCorrectionRel;
rec.minRsRelEigBefore = rsInfo.minRelEigBefore;
rec.rsCorrectionRel = rsInfo.correctionRel;

if covarianceGapRel > covGapTol
    rec.status = sprintf('covariance-gap-too-large-%.3e',covarianceGapRel);
    return;
end
if maxDesiredSignalGapRel > desiredGapTol
    rec.status = sprintf('desired-signal-gap-too-large-%.3e', ...
        maxDesiredSignalGapRel);
    return;
end
if maxSelfResidualLeakRel > leakTol
    rec.status = sprintf('self-residual-leak-too-large-%.3e', ...
        maxSelfResidualLeakRel);
    return;
end

rec.valid = true;
rec.status = 'covariance-preserving-recovery-ok';
end

function rec = make_empty_recovery(Nt,K)
rec.valid = false;
rec.status = 'not-run';
rec.failedUser = NaN;
rec.recoveryDenominator = NaN;
rec.recoveryDenominatorScale = NaN;
rec.W = [];
rec.Wc = zeros(Nt,K);
rec.Wr = [];
rec.RsRecovered = [];
rec.RxRecovered = [];
rec.RxSDR = [];
rec.residuals = zeros(Nt,Nt,K);
rec.power = NaN;
rec.powerSDR = NaN;
rec.covarianceGapRel = Inf;
rec.powerGapRel = Inf;
rec.maxDesiredSignalGapRel = Inf;
rec.maxSelfResidualLeakRel = Inf;
rec.minResidualRelEigBefore = NaN;
rec.maxResidualCorrectionRel = NaN;
rec.minRsRelEigBefore = NaN;
rec.rsCorrectionRel = NaN;
end

function [Apsd,info] = sanitize_psd(A,softTol,hardTol)
A = (A+A')/2;
info.valid = false;
info.minRelEigBefore = NaN;
info.correctionRel = Inf;
Apsd = A;
if any(~isfinite(A(:)))
    return;
end
[U,D] = eig(A);
d = real(diag(D));
if any(~isfinite(d))
    return;
end
scale = max(max(abs(d)),realmin);
minRel = min(d)/scale;
info.minRelEigBefore = minRel;
if minRel < -hardTol
    return;
end

% Only clip negative eigenvalues that are consistent with numerical error.
% Values between the soft and hard thresholds are still clipped, but the
% correction is recorded and can be rejected by the AO covariance-gap test.
dProj = d;
dProj(dProj < 0) = 0;
Apsd = U*diag(dProj)*U';
Apsd = (Apsd+Apsd')/2;
info.correctionRel = norm(Apsd-A,'fro')/max(norm(A,'fro'),realmin);
info.softProjectionUsed = minRel < -softTol;
info.valid = true;
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
