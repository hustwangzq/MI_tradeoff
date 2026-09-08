function [Jcrb, gradTheta, info, searchLogCRB] = SIM_crb_phase_value_gradient( ...
    params, ch, P, theta, Rx, alg)
%SIM_CRB_PHASE_VALUE_GRADIENT Exact CRB plus stabilized phase-search gradient.
%
% Jcrb and info always correspond to the exact unregularized physical CRB
% returned by SIM_crb_value. When useInternalCRBRegularization=true, the
% gradient is generated from the regularized CRB surrogate returned by
% SIM_crb_value_regularized. searchLogCRB is the scalar objective matched to
% that gradient. It is used only inside the phase-search model and is not
% stored as a reported CRB value.

N = params.N;
L = ch.L;

[Jcrb, info] = SIM_crb_value(ch, P, Rx, alg);
gradTheta = zeros(N,L);
gradLogTheta = zeros(N,L);
searchLogCRB = info.logCRB;
info.gradLogTheta = gradLogTheta;

if ~info.valid
    return;
end

Cx = info.Cx;
C0 = info.C0;
AinvGrad = info.Ainv;
BinvGrad = info.Binv;
aGrad = info.traceAinv;
bGrad = info.traceBinv;
isRRIS = isfield(ch,'type') && strcmpi(ch.type,'rRIS');
useNormalizedInverseGradient = isRRIS && ...
    isfield(info,'AinvNormalized') && ...
    ~isempty(info.AinvNormalized) && isfield(info,'BinvNormalized') && ...
    ~isempty(info.BinvNormalized);

useReg = logical(get_option(alg,'useInternalCRBRegularization',false));
regGradientActive = false;
regRelGradient = 0;
if useReg
    [~,regInfo] = SIM_crb_value_regularized(ch,P,Rx,alg);
    if regInfo.valid
        AinvGrad = regInfo.Ainv;
        BinvGrad = regInfo.Binv;
        aGrad = regInfo.traceAinv;
        bGrad = regInfo.traceBinv;
        searchLogCRB = regInfo.logCRB;
        regGradientActive = regInfo.regRel > 0;
        regRelGradient = regInfo.regRel;
    end
end

for ell = 1:L
    [UL, UR] = compute_UL_UR(theta, ch.Omega, ell);

    for n = 1:N
        En = zeros(N,N);
        En(n,n) = 1;
        dP = 1j*exp(1j*theta(n,ell)) * UL * En * UR;

        dA = dP*Cx*P' + P*Cx*dP';
        dB = dP*C0*P' + P*C0*dP';

        % If the optional trace-scaled loading is ever re-enabled, its
        % derivative must be included as well:
        %   d(A + eps*trace(A)/N*I)
        %       = dA + eps*trace(dA)/N*I.
        % Omitting this term makes the hand gradient inconsistent with the
        % regularized scalar objective near ill-conditioned points.
        if regGradientActive
            dA = dA + regRelGradient*real(trace(dA))/N*eye(N);
            dB = dB + regRelGradient*real(trace(dB))/N*eye(N);
        end

        da = -real(trace(AinvGrad*dA*AinvGrad));
        db = -real(trace(BinvGrad*dB*BinvGrad));
        if useNormalizedInverseGradient && ~regGradientActive
            daOverA = -real(trace(info.AinvNormalized*dA* ...
                info.AinvNormalized)) / max(info.maxEigA* ...
                info.traceAinvNormalized,realmin);
            dbOverB = -real(trace(info.BinvNormalized*dB* ...
                info.BinvNormalized)) / max(info.maxEigB* ...
                info.traceBinvNormalized,realmin);
            gradLog = daOverA + dbOverB;
        else
            gradLog = da/max(aGrad,realmin) + ...
                db/max(bGrad,realmin);
        end
        if isfinite(gradLog)
            gradLogTheta(n,ell) = real(gradLog);
        else
            gradLogTheta(n,ell) = 0;
        end

        % Compatibility output: original-objective scaling. The FP-SCA
        % phase routine uses gradLogTheta, not this scaled gradient.
        gradOriginal = Jcrb*gradLogTheta(n,ell);
        if isfinite(gradOriginal)
            gradTheta(n,ell) = real(gradOriginal);
        else
            gradTheta(n,ell) = 0;
        end
    end
end

info.gradLogTheta = gradLogTheta;
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
