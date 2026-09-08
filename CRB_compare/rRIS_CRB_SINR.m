% rRIS_CRB_SINR.m
% Run the stable covariance-domain CRB-SINR framework using the rRIS
% geometry/channel realization shared with MIc-Mis_revised/rRIS_main.m.
% The scenario flag is consumed and cleared by SIM_CRB_SINR.m.

setenv('CRB_SINR_SCENARIO','rRIS');
run(fullfile(fileparts(mfilename('fullpath')),'SIM_CRB_SINR.m'));
