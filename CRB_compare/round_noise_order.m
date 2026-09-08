function sigma2Clean = round_noise_order(sigma2Raw)
%ROUND_NOISE_ORDER Round a raw noise variance to a clean order of magnitude.
% Example: 1.8321e-10 -> 1e-10.  This keeps the intended SNR scale while
% avoiding long decimal noise values in the simulation output.

sigma2Raw = max(real(sigma2Raw), realmin);          % Avoid log10 of zero or negative numbers.
expo = floor(log10(sigma2Raw));                     % Keep the scientific-notation exponent.
sigma2Clean = 10^expo;                              % Return a clean one-digit order such as 1e-10.
end
