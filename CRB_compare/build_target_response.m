function [Hs, Rh] = build_target_response(risPos, targetPos, params)
%BUILD_TARGET_RESPONSE Point-target response and covariance prior.
%
% Paper correspondence:
%   H_s = sum_q alpha_q b_q b_q^H,
%   h_s = vec(H_s^H),
%   R_h = sum_q beta_q^2 (b_q^* \otimes b_q)(b_q^* \otimes b_q)^H.
%
% The target steering vector b_q is deterministic and unit-norm, as in the
% paper's UPA steering-vector model. The random target coefficient alpha_q
% carries the target reflection randomness and large-scale echo strength.
%
% Reproducibility rule:
%   If params.target.xi is provided, alpha_q = sqrt(beta_q^2)*xi_q is built
%   from these pre-generated CN(0,1) coefficients. This is the recommended
%   mode when comparing different L values, because target randomness is then
%   generated once outside the L-loop. If params.target.xi is absent, this
%   function falls back to deterministic per-target seeds for backward
%   compatibility.

N = params.N;                                      % Number of RIS/SIM elements per layer.
Qtar = size(targetPos,1);                          % Number of point targets.
lambda = params.lambda;                            % Carrier wavelength.
center = mean(risPos,1);                           % Center of the outermost RIS/SIM layer.
Hs = zeros(N,N);                                   % Instantaneous target response matrix H_s.
Rh = zeros(N^2,N^2);                               % Prior covariance R_{h_s}.
sigmaAlpha2 = params.target.sigmaAlpha2;           % Target reflection variance before pathloss scaling.

usePreGeneratedXi = isfield(params,'target') ...    % Whether common target coefficients are supplied.
    && isfield(params.target,'xi') ...
    && numel(params.target.xi) >= Qtar;

for q = 1:Qtar                                     % Add each point-target contribution.
    bq = farfield_steering_vector(risPos, targetPos(q,:), lambda, true); % Unit-norm target steering vector b_q.
    d0 = norm(targetPos(q,:) - center);             % Outermost-layer-to-target distance.
    pl = pathloss_linear(d0, params);               % One-way large-scale pathloss for target direction.
    beta2 = (pl^2) * sigmaAlpha2;                   % Variance beta_q^2 of alpha_q, using two-way echo scaling.

    if usePreGeneratedXi                            % Preferred: reuse common target realization.
        xi_q = params.target.xi(q);                 % Fixed CN(0,1) coefficient generated outside this function.
    else                                           % Backward-compatible fallback.
        oldStream = rng;                            % Store random state before target coefficient draw.
        rng(params.seed.targetBase + q);            % Deterministic target coefficient for this q.
        xi_q = (randn + 1j*randn)/sqrt(2);          % CN(0,1) coefficient.
        rng(oldStream);                             % Restore the previous random state.
    end

    alpha = sqrt(beta2) * xi_q;                     % Instantaneous complex coefficient alpha_q.
    Hs = Hs + alpha * (bq*bq');                     % H_s contribution alpha_q b_q b_q^H.
    vq = kron(conj(bq), bq);                        % vec((b_q b_q^H)^H) basis = b_q^* \otimes b_q.
    Rh = Rh + beta2 * (vq*vq');                     % Add beta_q^2 v_q v_q^H to R_h.
end

loadScale = params.target.rhDiagLoad * max(real(trace(Rh))/max(N^2,1), 1e-30); % Tiny background loading scale.
Rh = Rh + loadScale*eye(N^2);                       % Regularize low-rank point-target covariance slightly.
Rh = nearest_hermitian_pd(Rh, 1e-14);               % Final numerical Hermitian positive-definite cleanup.
end
