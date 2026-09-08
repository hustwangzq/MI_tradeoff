function [sigma_c2, sigma_s2] = set_noise_powers(params, ch, P, W, S)
%SET_NOISE_POWERS Set clean communication and sensing noise variances.
% The raw noise powers are computed from an initial received power and target
% SNR. They are optionally rounded to their order of magnitude, e.g.,
% 1.83e-10 -> 1e-10, so simulation settings are easy to report.

K = params.K;                                      % Number of communication users.
pc = 0;                                            % Accumulator for communication received stream power.
for k = 1:K                                        % Average over users.
    yrow = ch.hUsers{k}' * P * ch.G * W;            % Effective stream gains at CU k.
    pc = pc + sum(abs(yrow).^2);                    % Total stream power at CU k before noise.
end
pc = pc / K;                                       % Average communication signal power.
sigma_c2_raw = max(pc/10^(params.SNRc_dB/10), 1e-30); % Raw communication noise variance.

% X = W*S;                                           % Transmit signal matrix X = W S.
% Y0 = ch.G' * P' * ch.Hs * P * ch.G * X;             % Noise-free sensing echo at the BS.
% ps = norm(Y0,'fro')^2 / max(numel(Y0),1);           % Average sensing echo power per received sample.
% sigma_s2_raw = max(ps/10^(params.SNRs_dB/10), 1e-30); % Raw sensing noise variance.

B = sensing_B_matrix(ch, P, W, S);                  % Sensing observation matrix B_s.
ps = real(trace(B*ch.Rh*B')) / max(size(B,1),1);    % Average statistical sensing echo power.
sigma_s2_raw = max(ps/10^(params.SNRs_dB/10), 1e-30); % Raw sensing noise variance.


if isfield(params,'noise') && params.noise.roundToOrder % If clean powers are requested,
    sigma_c2 = round_noise_order(sigma_c2_raw);     % round communication noise to nearest order.
    sigma_s2 = round_noise_order(sigma_s2_raw);     % round sensing noise to nearest order.
else                                                % Otherwise,
    sigma_c2 = sigma_c2_raw;                        % keep exact communication noise variance.
    sigma_s2 = sigma_s2_raw;                        % keep exact sensing noise variance.
end
end
