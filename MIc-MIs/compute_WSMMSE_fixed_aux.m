function J = compute_WSMMSE_fixed_aux(params, ch, P, W, S, aux)
%COMPUTE_WSMMSE_FIXED_AUX Positive fixed-auxiliary WSMMSE objective.
%
% Paper correspondence:
%   This function evaluates the objective of the phase subproblem (31):
%       J = w_c sum_k epsilon_k E_{c,k} + (w_s/N) tr(V_s E_s),
%   with f_k, F_s, epsilon_k, and V_s fixed.
%
% Use in code:
%   It is used by the Armijo backtracking line search. Therefore E_{c,k}
%   must be evaluated by the full formula (15), not by the simplified
%   optimal-MSE expression 1-|a_k|^2/p_k, because f_k is fixed during the
%   line search.

K = params.K;                                      % Number of CUs.
N = params.N;                                      % Number of RIS/SIM elements per layer.
wc = params.wc;                                    % Communication weight w_c.
ws = params.ws;                                    % Sensing weight w_s.

Gc = zeros(params.Nt,K);                            % Equivalent communication channels g_k.
for k = 1:K                                        % Loop over CUs.
    Gc(:,k) = ch.G' * P' * ch.hUsers{k};            % g_k = G_st^H P^H h_k.
end
Auser = Gc' * W;                                    % a_{k,m}=g_k^H w_m for all streams.

Jc = 0;                                             % Communication contribution to J.
for k = 1:K                                        % Loop over CUs.
    p_k = sum(abs(Auser(k,:)).^2) + ch.sigma_c2;    % p_k = ||g_k^H W||^2 + sigma_c^2.
    Ec_k = abs(aux.f(k))^2*p_k ...                  % |f_k|^2 p_k term in (15).
         - 2*real(aux.f(k)*Auser(k,k)) + 1;         % -2 Re{f_k g_k^H w_k}+1 term in (15).
    Jc = Jc + params.userWeight(k)*aux.epsilon(k)*real(Ec_k); % Sum epsilon_k E_{c,k}.
end

B = sensing_B_matrix(ch, P, W, S);                  % Current sensing matrix B_s.
Fs = aux.Fs;                                        % Fixed sensing filter F_s.
Vs = aux.Vs;                                        % Fixed sensing weight V_s.
Rh = ch.Rh;                                         % Target-response covariance R_h.
Es = (Fs*B-eye(N^2))*Rh*(Fs*B-eye(N^2))' ...        % E_s with fixed F_s.
   + ch.sigma_s2*(Fs*Fs');                          % Noise term sigma_s^2 F_s F_s^H.
Js = real(trace(Vs*Es));                            % Sensing contribution tr(V_s E_s).

J = wc*Jc + (ws/N)*Js;                              % Weighted positive WSMMSE objective.
end
