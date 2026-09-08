function [gradGTheta, J0] = phase_gradient_WSMMSE(params, ch, P, theta, W, S, aux)
%PHASE_GRADIENT_WSMMSE Phase gradient for the fixed-auxiliary WSMMSE subproblem.
%
% Paper/code correspondence:
%   Stage V minimizes the positive fixed-auxiliary WSMMSE phase objective
%
%       J_phi(theta) =
%           w_c sum_k epsilon_k E_c,k
%           + (w_s/N) tr(V_s E_s),
%
%   where f_k, F_s, epsilon_k, V_s, and W are fixed during the phase update.
%
%   Appendix A defines g = -J_phi for the phase variables. Therefore,
%
%       grad_g = d g / d theta = - dJ_phi / dtheta.
%
% Important implementation note:
%   The Armijo line search evaluates J_phi by compute_WSMMSE_fixed_aux().
%   Therefore, the sensing gradient here is computed as the direct derivative
%   of tr(V_s E_s) with fixed F_s and V_s, instead of using the Appendix-A
%   E_s V_s E_s equivalent-gradient expression. This guarantees that the
%   search direction matches the Armijo objective exactly.
%
% Output:
%   gradGTheta(n,ell) = nabla_{theta_{ell,n}} g.
%   The phase update theta <- theta + alpha*gradGTheta is equivalent to
%   gradient descent on the positive objective J_phi.

N  = params.N;                                      % Number of elements per layer.
K  = params.K;                                      % Number of communication users.
M  = params.Ns;                                     % Total streams M=K+N_t.
L  = ch.L;                                          % Number of RIS/SIM layers.
wc = params.wc;                                     % Communication weight w_c.
ws = params.ws;                                     % Sensing weight w_s.

if isfield(params,'userWeight')                     % Optional per-user weights.
    userWeight = params.userWeight(:);              % K-by-1 user weights.
else                                                % If not provided,
    userWeight = ones(K,1);                         % use equal weights.
end

% Current stream gains a_{k,m}=h_k^H P G W.
Auser = zeros(K,M);                                  % Row k stores all stream gains of CU k.
for k = 1:K                                         % Loop over all CUs.
    Auser(k,:) = ch.hUsers{k}' * P * ch.G * W;       % h_k^H P G_st W.
end

% Quantities for the direct fixed-auxiliary sensing derivative.
X = W*S;                                            % Transmit signal matrix X = W S.
B = sensing_B_matrix(ch, P, W, S);                  % Current sensing matrix B_s.
Cerr = aux.Fs*B - eye(N^2);                         % F_s B_s - I in E_s.
Rh = ch.Rh;                                         % Target covariance R_h.
Fs = aux.Fs;                                        % Fixed sensing receive filter.
Vs = aux.Vs;                                        % Fixed sensing WMMSE weight.

% Fixed factors in B_s = kron((P G)^T, X^H G^H P^H).
Cbase = (P*ch.G).';                                 % (P G_st)^T.
Mbase = X' * ch.G' * P';                            % X^H G_st^H P^H.

J0 = compute_WSMMSE_fixed_aux(params, ch, P, W, S, aux); % Positive objective J_phi for Armijo.
gradGTheta = zeros(N,L);                             % Allocate nabla_theta g.

for ell = 1:L                                       % Loop over all layers.
    [UL, UR] = compute_UL_UR(theta, ch.Omega, ell);  % P = UL * Phi_ell * UR.

    for n = 1:N                                     % Loop over all elements in layer ell.
        En = zeros(N,N);                            % Selection matrix e_n e_n^T.
        En(n,n) = 1;                                % Select element n.

        dPhi = 1j*exp(1j*theta(n,ell))*En;           % d Phi_ell / d theta_{ell,n}.
        dP = UL * dPhi * UR;                        % dP / d theta_{ell,n}.

        %% ---------------- Communication derivative of J_phi ----------------
        dJ_c = 0;                                   % Communication contribution d/dtheta sum eps E_c.
        for k = 1:K                                 % Sum over all users.
            da = ch.hUsers{k}' * dP * ch.G * W;      % d a_{k,m} for all streams m.
            dp = 2*real(sum(conj(Auser(k,:)).*da)); % d p_k = 2 Re{sum_m a_{k,m}^* d a_{k,m}}.

            dEc = abs(aux.f(k))^2 * dp ...          % |f_k|^2 d p_k.
                - 2*real(aux.f(k) * da(k));         % -2 Re{f_k d a_{k,k}}.

            dJ_c = dJ_c + userWeight(k)*aux.epsilon(k)*dEc; % Weighted dE_c,k.
        end

        %% ---------------- Sensing derivative of J_phi ----------------------
        % Direct derivative of B_s:
        %
        %   B_s = kron((P G)^T, X^H G^H P^H).
        %
        % Since X, F_s, V_s, and R_h are fixed in Stage V,
        %
        %   dB_s = kron((dP G)^T, X^H G^H P^H)
        %        + kron((P G)^T, X^H G^H dP^H).
        %
        dB = kron((dP*ch.G).', Mbase) ...           % d(PG)^T \otimes X^H G^H P^H.
           + kron(Cbase, X' * ch.G' * dP');         % (PG)^T \otimes X^H G^H dP^H.

        % E_s = (F_s B_s - I) R_h (F_s B_s - I)^H + sigma_s^2 F_s F_s^H.
        %
        % With fixed F_s and V_s,
        %
        % d tr(V_s E_s)
        %   = 2 Re{ tr( V_s F_s dB_s R_h (F_s B_s - I)^H ) }.
        %
        dJ_s = 2*real(trace(Vs * Fs * dB * Rh * Cerr')); % Direct fixed-aux sensing gradient.

        %% ---------------- Total gradient and negative-WSMMSE gradient ------
        dJ = wc*dJ_c + (ws/N)*dJ_s;                 % dJ_phi / d theta_{ell,n}.
        gradGTheta(n,ell) = -dJ;                    % nabla_theta g = - dJ_phi/dtheta.
    end
end
end