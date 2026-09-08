function [cScaled, gradScaled, details] = SIM_sinr_phase_constraints( ...
    params, ch, P, theta, Qcomm, Rs, gammaTarget, scaleIn, needGradient)
%SINR_PHASE_CONSTRAINTS Evaluate phase-dependent SINR inequalities.
%
% The raw inequality for user k is
%   c_k = gamma_k*(interference_k + noise) - desired_k <= 0.
%
% To improve numerical scaling, c_k and its gradient
% are divided by a fixed positive scale. If scaleIn is empty, the scale is
% constructed from the current point and returned in details.scale. During
% one phase block, the caller should reuse this fixed scale.

if nargin < 8
    scaleIn = [];
end
if nargin < 9 || isempty(needGradient)
    needGradient = true;
end

K = params.K;
N = params.N;
L = ch.L;
gammaTarget = gammaTarget(:);
Rx = Rs;
for j = 1:K, Rx = Rx + Qcomm(:,:,j); end
Rx = (Rx+Rx')/2;

raw = zeros(K,1);
desired = zeros(K,1);
interferenceNoise = zeros(K,1);
for k = 1:K
    gk = ch.G' * P' * ch.hUsers{k};
    totalPower = real(gk'*Rx*gk) + ch.sigma_c2;
    desired(k) = real(gk'*Qcomm(:,:,k)*gk);
    interferenceNoise(k) = totalPower - desired(k);
    raw(k) = gammaTarget(k)*interferenceNoise(k) - desired(k);
end

if isempty(scaleIn)
    scale = abs(gammaTarget.*interferenceNoise) + abs(desired);
    scale = max(scale, gammaTarget*ch.sigma_c2);
    scale = max(scale, 1e-30);
else
    scale = scaleIn(:);
end

cScaled = raw ./ scale;
gradScaled = [];

if needGradient
    gradScaled = zeros(N,L,K);

    for ell = 1:L
        [UL, UR] = compute_UL_UR(theta, ch.Omega, ell);

        for n = 1:N
            En = zeros(N,N);
            En(n,n) = 1;
            dP = 1j*exp(1j*theta(n,ell)) * UL * En * UR;

            for k = 1:K
                hk = ch.hUsers{k};
                gk = ch.G' * P' * hk;
                dg = ch.G' * dP' * hk;
                dDesired = 2*real(dg'*Qcomm(:,:,k)*gk);
                dTotal = 2*real(dg'*Rx*gk);
                dInterference = dTotal - dDesired;
                dRaw = gammaTarget(k)*dInterference - dDesired;
                gradScaled(n,ell,k) = dRaw/scale(k);
            end
        end
    end
end

details.raw = raw;
details.scale = scale;
details.desired = desired;
details.interferenceNoise = interferenceNoise;
details.sinr = desired ./ max(interferenceNoise,1e-30);
details.maxScaledViolation = max([0; cScaled]);
end
