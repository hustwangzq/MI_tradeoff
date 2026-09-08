function B = sensing_B_matrix(ch, P, W, S)
%SENSING_B_MATRIX Build B_s for y_s = vec(Y_s^H) = B_s vec(H_s^H)+z_s.
%
% Sensing model:
%   Y_s = G_st^H P^H H_s P G_st X + Z_s,  X = W S.
% After taking the Hermitian transpose and vectorizing,
%   vec(Y_s^H) = [(P G_st)^T \otimes (X^H G_st^H P^H)] vec(H_s^H) + vec(Z_s^H).

X = W*S;                                           % Transmit signal matrix X = W S, size N_t-by-T.
C = (P*ch.G).';                                    % First Kronecker factor (P G_st)^T, size N_t-by-N.
M = X' * ch.G' * P';                               % Second factor X^H G_st^H P^H, size T-by-N.
B = kron(C, M);                                    % B_s, size (T N_t)-by-N^2.
end
