function [sigma_c2_ref, sigma_s2_ref, noiseInfo] = estimate_reference_noise_average(params, chRef, seedBase, numTrials)
%ESTIMATE_REFERENCE_NOISE_AVERAGE Estimate fixed reference noise powers.
%
% This function estimates one common noise floor from a reference channel by
% averaging over multiple random initial W, S, and theta. The returned noise
% powers are used for all compared L values in the same scenario.
%
% Important:
%   roundToOrder is not used inside each trial. Averaging should be performed
%   on raw values. Optional lower bounds are applied after averaging.

if nargin < 4 || isempty(numTrials)
    numTrials = 20;
end

paramsNoise = params;                              % Local copy for noise estimation.
paramsNoise.noise.roundToOrder = true;            % Average raw noise values, not rounded values.

sigma_c2_list = zeros(numTrials,1);                % Communication noise samples.
sigma_s2_list = zeros(numTrials,1);                % Sensing noise samples.

for tt = 1:numTrials
    initTmp = initialize_solver_state(paramsNoise, chRef, seedBase + tt);
    [Ptmp, ~] = build_P(initTmp.theta0, chRef.Omega);

    [sigma_c2_list(tt), sigma_s2_list(tt)] = ...
        set_noise_powers(paramsNoise, chRef, Ptmp, initTmp.W0, initTmp.S);
end

sigma_c2_raw_avg = mean(sigma_c2_list);            % Average raw communication noise.
sigma_s2_raw_avg = mean(sigma_s2_list);            % Average raw sensing noise.

if isfield(params,'noise') && isfield(params.noise,'minSigmaC2')
    minSigmaC2 = params.noise.minSigmaC2;
else
    minSigmaC2 = 1e-30;
end

if isfield(params,'noise') && isfield(params.noise,'minSigmaS2')
    minSigmaS2 = params.noise.minSigmaS2;
else
    minSigmaS2 = 1e-30;
end

sigma_c2_ref = max(sigma_c2_raw_avg, minSigmaC2);  % Apply communication noise floor.
sigma_s2_ref = max(sigma_s2_raw_avg, minSigmaS2);  % Apply sensing noise floor.

if isfield(params,'noise') && params.noise.roundToOrder
    sigma_c2_ref = round_noise_order(sigma_c2_ref);
    sigma_s2_ref = round_noise_order(sigma_s2_ref);
end

noiseInfo.sigma_c2_list = sigma_c2_list;
noiseInfo.sigma_s2_list = sigma_s2_list;
noiseInfo.sigma_c2_raw_avg = sigma_c2_raw_avg;
noiseInfo.sigma_s2_raw_avg = sigma_s2_raw_avg;
noiseInfo.sigma_c2_ref = sigma_c2_ref;
noiseInfo.sigma_s2_ref = sigma_s2_ref;
noiseInfo.minSigmaC2 = minSigmaC2;
noiseInfo.minSigmaS2 = minSigmaS2;
noiseInfo.numAvg = numTrials;       
end