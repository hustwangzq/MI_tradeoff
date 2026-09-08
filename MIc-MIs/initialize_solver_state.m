function state = initialize_solver_state(params, ch, seed)
%INITIALIZE_SOLVER_STATE Generate a common feasible initialization.
% The same initialization is reused across the weight sweep so that the MI
% tradeoff curve is not distorted by different random starting points.

if nargin >= 3 && ~isempty(seed)                    % If a seed is supplied,
    rng(seed);                                      % fix the random generator for reproducibility.
end

Nt = params.Nt;                                     % Number of BS antennas.
Ns = params.Ns;                                     % Total number of streams, Ns=K+Nt.
N  = params.N;                                      % Number of elements per layer.
L  = ch.L;                                          % Number of RIS/SIM layers.
T  = params.T;                                      % Number of sensing snapshots.

state.S = (randn(Ns,T) + 1j*randn(Ns,T))/sqrt(2);   % Stream matrix S with unit-power Gaussian entries.
state.W0 = (randn(Nt,Ns)+1j*randn(Nt,Ns))/sqrt(2);  % Random initial transmit beamforming matrix W.
state.W0 = sqrt(params.P0) * state.W0 / norm(state.W0,'fro'); % Normalize W to satisfy tr(W W^H)=P0.
state.theta0 = 2*pi*rand(N,L);                      % Random initial phases theta_{ell,n} in [0,2pi).
% [P0, ~] = build_P(state.theta0, ch.Omega);           % Build initial end-to-end matrix P from theta0.
% [state.sigma_c2, state.sigma_s2] = set_noise_powers(params, ch, P0, state.W0, state.S); % Set clean noise powers.
end
