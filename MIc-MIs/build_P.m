function [P, Phi] = build_P(theta, Omega)
%BUILD_P Construct P = Phi_L Omega_L ... Phi_2 Omega_2 Phi_1.
% WZQ SHUAI
[N,L] = size(theta);
Phi = cell(L,1);
for ell = 1:L
    Phi{ell} = diag(exp(1j*theta(:,ell)));
end
P = Phi{1};
for ell = 2:L
    P = Phi{ell} * Omega{ell} * P;
end
end
