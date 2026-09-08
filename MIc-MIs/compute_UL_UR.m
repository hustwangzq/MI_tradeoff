function [UL, UR] = compute_UL_UR(theta, Omega, ell)
%COMPUTE_UL_UR Decompose P as P = UL * Phi_ell * UR.
%
% Paper correspondence:
%   \dot P_{ell,n} = j exp(j theta_{ell,n}) UL e_n e_n^T UR.
%
% The implemented end-to-end matrix is
%   P = Phi_L Omega_L Phi_{L-1} ... Omega_2 Phi_1.
% Therefore,
%   UR = Omega_ell Phi_{ell-1} ... Omega_2 Phi_1, for ell>1,
%   UL = Phi_L Omega_L ... Phi_{ell+1} Omega_{ell+1}, for ell<L.

[N,L] = size(theta);                                % N elements and L layers.
Phi = cell(L,1);                                    % Phase-shift matrices Phi_ell.
for m = 1:L                                        % Build all phase-shift matrices.
    Phi{m} = diag(exp(1j*theta(:,m)));              % Phi_m = diag(exp(j theta_{m,n})).
end

if ell == 1                                        % If Phi_ell is the first layer,
    UR = eye(N);                                    % there is no right-side product.
else                                                % Otherwise build right-side product.
    UR = Phi{1};                                    % Start from Phi_1.
    for m = 2:ell-1                                % Add Omega_m and Phi_m up to ell-1.
        UR = Phi{m} * Omega{m} * UR;                % UR = Phi_m Omega_m ... Phi_1.
    end
    UR = Omega{ell} * UR;                           % Include Omega_ell immediately right of Phi_ell.
end

UL = eye(N);                                        % Initialize left-side product.
for m = L:-1:(ell+1)                               % Multiply layers from L down to ell+1.
    UL = UL * Phi{m} * Omega{m};                    % UL = Phi_L Omega_L ... Phi_{ell+1} Omega_{ell+1}.
end
end
