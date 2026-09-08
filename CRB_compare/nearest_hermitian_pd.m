function A = nearest_hermitian_pd(A, relJitter)
%NEAREST_HERMITIAN_PD Small numerical Hermitian positive-definite regularization.
%
% This routine is NOT part of the theoretical algorithm. It only prevents
% numerical failures when inverting matrices that should be Hermitian positive
% definite in theory but may have tiny negative eigenvalues due to roundoff.
% The diagonal loading is relative to the matrix norm and is intentionally tiny.

if nargin < 2 || isempty(relJitter)                 % If no relative jitter is provided,
    relJitter = 1e-12;                              % use a tiny relative diagonal loading.
end

A = (A + A')/2;                                     % Enforce Hermitian symmetry first.
scaleA = max(1, norm(A,'fro'));                     % Reference scale of the matrix.
jitter = relJitter * scaleA;                        % Absolute jitter adapted to the matrix magnitude.

[R,p] = chol(A + jitter*eye(size(A)));              % Check whether a tiny loading is sufficient.
if p == 0                                          % If Cholesky succeeds,
    A = A + jitter*eye(size(A));                    % keep the minimally loaded Hermitian matrix.
    return;                                         % Finish without expensive eigenvalue correction.
end

eigMin = min(real(eig(A)));                         % Smallest eigenvalue used only if Cholesky failed.
if eigMin < jitter                                  % If the matrix is not numerically positive definite,
    A = A + (jitter - eigMin)*eye(size(A));          % shift all eigenvalues just above the jitter level.
else                                                % Otherwise,
    A = A + jitter*eye(size(A));                    % add a tiny diagonal loading for stable inversion.
end
end
