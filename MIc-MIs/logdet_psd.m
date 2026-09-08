function val = logdet_psd(A)
%LOGDET_PSD Stable log-det for Hermitian positive semidefinite matrices.

A = (A+A')/2;
jitter = 1e-12;
[R,p] = chol(A + jitter*eye(size(A)));
if p == 0
    val = 2*sum(log(real(diag(R))));
else
    eigv = max(real(eig(A)), 1e-12);
    val = sum(log(eigv));
end
end
