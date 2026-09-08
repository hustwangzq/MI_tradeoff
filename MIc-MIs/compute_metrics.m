function metrics = compute_metrics(params, ch, P, W, S)
%COMPUTE_METRICS Communication MI, sensing MI, sum rate, and NMMSE.
% Communication MI is sum_k log(1+SINR_k). Sensing MI is computed from
% I_s = log det(I + sigma_s^{-2} B_s R_h B_s^H).

K = params.K;                                      % Number of communication users.
N = params.N;                                      % Number of RIS/SIM elements per layer.
MIc = 0;                                           % Communication MI accumulator.
rate = 0;                                          % Sum-rate accumulator in bit/s/Hz.
for k = 1:K                                        % Evaluate every user.
    gk = ch.G' * P' * ch.hUsers{k};                % Equivalent channel g_k.
    a = gk' * W;                                   % Stream gains g_k^H W.
    pTotal = sum(abs(a).^2) + ch.sigma_c2;         % Desired + interference + sensing streams + noise.
    desired = abs(a(k))^2;                         % Desired stream power.
    interf = pTotal - desired;                     % Interference-plus-noise power.
    gamma = desired/max(interf,1e-16);             % SINR of user k.
    MIc = MIc + log(1+gamma);                      % Natural-log communication MI contribution.
    rate = rate + log2(1+gamma);                   % Bit/s/Hz rate contribution.
end

B = sensing_B_matrix(ch, P, W, S);                  % Sensing observation matrix B_s.
Rh = ch.Rh;                                         % Target-response covariance R_h.
MIs = real(logdet_psd(eye(size(B,1)) + (1/ch.sigma_s2)*(B*Rh*B'))); % Sensing MI.

RhReg = nearest_hermitian_pd(Rh, 1e-12);            % Regularized covariance for posterior-MMSE diagnostic.
PostInfo = (RhReg \ eye(N^2)) + (1/ch.sigma_s2)*(B'*B); % R_h^{-1}+sigma_s^{-2}B_s^H B_s.
PostInfo = nearest_hermitian_pd(PostInfo, 1e-12);   % Ensure numerical positive definiteness.
Epost = PostInfo \ eye(N^2);                        % Posterior covariance used only for NMMSE monitoring.
Epost = (Epost+Epost')/2;                           % Remove tiny asymmetry.
nmmse = real(trace(Epost))/max(real(trace(Rh)),1e-15); % Normalized sensing MMSE.

metrics.MIc = MIc;                                  % Communication MI.
metrics.MIs = MIs;                                  % Sensing MI.
metrics.rate = rate;                                % Communication sum rate in bit/s/Hz.
metrics.nmmse = nmmse;                              % Normalized sensing MMSE.
end
