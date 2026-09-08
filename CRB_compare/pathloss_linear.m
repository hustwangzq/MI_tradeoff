function pl = pathloss_linear(dist, params)
%PATHLOSS_LINEAR Large-scale power gain with 1 m reference distance.

dist = max(dist, 1e-12);
PLdB = params.pathlossRef_dB - 10*params.pathlossExp*log10(dist);
pl = 10^(PLdB/10);
end
