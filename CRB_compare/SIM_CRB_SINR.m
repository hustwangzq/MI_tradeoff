clear; close all; clc;

% SIM_CRB_SINR.m
% SIM-only CRB-SINR tradeoff simulation.
%
% Main revisions in this separated version:
%   1) every SIM optimization function uses the SIM_ prefix;
%   2) the first formal point of each L is selected by CRB-oriented
%      multi-initialization, not by the maximum-SINR range criterion;
%   3) every retained CRB value is stored together with the W, theta, P, Rx,
%      and all metrics that produced it;
%   4) the 1-dB coarse grid is retained, while abnormal intervals receive
%      complete linear-SINR midpoint insertion with a bounded depth;
%   5) remote random/structured branch switching is delayed until local
%      refinement is exhausted;
%   6) the phase block uses feasible monotone proximal SCA: every accepted
%      phase satisfies the exact SINR test and does not increase exact CRB;
%   7) selected L-dependent SINR windows are monitored using four
%      linear-SINR subintervals and a larger AO/phase iteration budget;
%   8) every formal SINR point propagates two independent states: the
%      reported best-CRB solution and a communication-reserve carrier;
%      carriers are saved separately and are never plotted as Pareto points;
%   9) backward re-optimization is disabled in this diagnostic version so
%      the plotted curve is a pure forward dual-state continuation path;
%  10) two figures are drawn: exact unregularized CRB and a transient
%      regularized-CRB diagnostic evaluated from the same retained solutions.
%      Regularized CRB values are not written into the MAT result records.

scriptVersion = 'SIM-FP-SCA-persistent-bestCRB-commCarrier-v12-diagnostics-20260825';

thisFile = mfilename('fullpath');
if isempty(thisFile)
    packageDir = pwd;
else
    packageDir = fileparts(thisFile);
end
scenarioCRB=getenv('CRB_SINR_SCENARIO');
setenv('CRB_SINR_SCENARIO','');
if isempty(scenarioCRB), scenarioCRB='SIM'; end
isRRIS=strcmpi(scenarioCRB,'rRIS');
if isRRIS
    scenarioCRB='rRIS';
    scriptVersion=['rRIS-' scriptVersion];
elseif ~strcmpi(scenarioCRB,'SIM')
    error('Unsupported CRB-SINR scenario: %s',scenarioCRB);
else
    scenarioCRB='SIM';
end
if strcmp(getenv('SIM_CRB_DIAGNOSTIC_SELFTEST'),'1')
    run_transition_diagnostic_selftest(packageDir,scriptVersion);
    return;
end
addpath(packageDir,'-begin');
rehash path;

if exist('cvx_begin','file') ~= 2
    error(['CVX is required. Install CVX, run cvx_setup, and rerun ', ...
        'SIM_CRB_SINR.m.']);
end

requiredSIMFunctions = { ...
    'SIM_build_common_links', ...
    'SIM_build_channels', ...
    'SIM_build_cross_layer_inits', ...
    'SIM_search_max_common_sinr_AO', ...
    'SIM_run_CRB_SINR_AO_solver', ...
    'SIM_recover_CRB_beamformers', ...
    'SIM_build_comm_carrier', ...
    'SIM_solve_comm_margin_SDR', ...
    'SIM_update_all_phase_margin', ...
    'SIM_print_transition_diagnostics', ...
    'SIM_solve_CRB_SINR_SDR', ...
    'SIM_solve_CRB_SINR_SDR_robust', ...
    'SIM_build_adaptive_CRB_regularization_plan', ...
    'SIM_estimate_pure_comm_sinr_boundary', ...
    'SIM_max_common_sinr_fixed_phase_SDR', ...
    'SIM_build_boundary_aware_gamma_grid', ...
    'SIM_rescue_high_sinr_target', ...
    'SIM_update_outer_phase_margin', ...
    'SIM_update_CRB_phase_FP_SCA', ...
    'SIM_solve_phase_FP_SCA_subproblem', ...
    'SIM_crb_phase_value_gradient', ...
    'SIM_sinr_phase_constraints', ...
    'SIM_evaluate_CRB_SINR_metrics', ...
    'SIM_crb_value', ...
    'SIM_crb_value_regularized'};
if isRRIS
    requiredSIMFunctions=[requiredSIMFunctions, ...
        {'rRIS_build_common_bank','rRIS_build_channels'}];
end
for ii = 1:numel(requiredSIMFunctions)
    activeFile = which(requiredSIMFunctions{ii});
    if isempty(activeFile)
        error('Missing required SIM function: %s.m',requiredSIMFunctions{ii});
    end
    if ~strcmpi(fileparts(activeFile),packageDir)
        error(['Path conflict for %s. Active file:\n%s\n', ...
            'Expected folder:\n%s'],requiredSIMFunctions{ii}, ...
            activeFile,packageDir);
    end
end

fprintf('[%s] Script version: %s\n',scenarioCRB,scriptVersion);
fprintf('[%s] Package folder: %s\n',scenarioCRB,packageDir);
fprintf('[CRB-SINR] Scenario: %s\n',scenarioCRB);

params0 = default_params();
if params0.Nt < params0.N
    error('The full response-matrix CRB requires Nt >= N.');
end
rng(params0.seed.global);

L_list = 1:3;
modeSIM = 1;
numNoiseAvg = 30;

% Conservative formal low-SINR targets.  SIM keeps its validated grid.
% For rRIS, a fixed-phase pure-communication SDR sampled on the locked MI
% phases gives representative common-SINR ceilings of about 30.36, 14.17,
% and -0.02 dB for L=1,2,3.  The rRIS curve therefore starts well below
% those ceilings at 0,-5,-10 dB; feasibility recovery can still move lower.
gammaTradeStartSIM = [0.5,0.1,0.05];
if isRRIS
    gammaTradeStartSIM = 10.^([0,-5,-10]/10);
end
gammaSeedSIM = gammaTradeStartSIM;

% ------------------------- Range-search options -------------------------
searchOptsMultiplicative.searchMode = 'multiplicative';
searchOptsMultiplicative.growthFactor = 1.3;
searchOptsMultiplicative.maxOuterGrowth = 40;
searchOptsMultiplicative.maxOuterConfirm = 100;
searchOptsMultiplicative.maxOuterRefine = 60;
searchOptsMultiplicative.maxExpand = 100;
searchOptsMultiplicative.maxBridge = 3;
searchOptsMultiplicative.maxRefine = 8;
searchOptsMultiplicative.rangeRelTol = 0.01;
searchOptsMultiplicative.enableRecovery = false;
searchOptsMultiplicative.enableFullBudgetBoundary = true;
searchOptsMultiplicative.maxBoundaryTrials = 5;
searchOptsMultiplicative.boundaryStepRel = 0.01;
searchOptsMultiplicative.boundaryStepMinRel = 1e-4;
searchOptsMultiplicative.boundaryStepMinAbs = 1e-4;
searchOptsMultiplicative.boundaryStepShrink = 0.5;

searchOptsAdditive.searchMode = 'additive';
searchOptsAdditive.maxOuterGrowth = 50;
searchOptsAdditive.maxOuterConfirm = 120;
searchOptsAdditive.maxOuterRefine = 80;
searchOptsAdditive.maxExpand = 120;
searchOptsAdditive.maxBridge = 0;
searchOptsAdditive.maxRefine = 10;
searchOptsAdditive.rangeRelTol = 0.005;
searchOptsAdditive.enableRecovery = true;
searchOptsAdditive.recoveryFactor = 0.90;
searchOptsAdditive.maxRecovery = 35;
searchOptsAdditive.stepShrink = 0.5;
searchOptsAdditive.maxStepHalving = 10;
searchOptsAdditive.enableFullBudgetBoundary = true;
searchOptsAdditive.maxBoundaryTrials = 5;
searchOptsAdditive.boundaryStepRel = 0.01;
searchOptsAdditive.boundaryStepMinRel = 1e-4;
searchOptsAdditive.boundaryStepMinAbs = 1e-4;
searchOptsAdditive.boundaryStepShrink = 0.5;

additiveStepRel = 0.02;
additiveStepAbsMin = 0.05;
additiveStepMinRel = 1e-3;
additiveStepMinAbs = 1e-3;

% -------------------------- Formal-curve options ------------------------
% Keep the original 1-dB coarse grid, but refine only abnormal intervals.
formalStepDB = 1.0;
backwardPolishRelTol = 1e-10; % legacy compatibility only

% No direct backward envelope is used. Instead, a higher-SINR complete
% solution may be used once as a genuine AO initialization for the adjacent
% lower target when the raw forward curve bends left by more than 0.10 dB.
enableBackwardPolish = false;
enableMonotoneEnvelope = false;
monotoneEnvelopeRelTol = 1e-10;
enableBackwardReoptimization = false;
backwardRepairTolDB = 0.10; % compatibility/diagnostic field only
maxBackwardRepairPass = 0;

% ---------------------- Dual-state continuation -------------------------
% Every retained Pareto point remains the report solution that minimizes
% exact physical CRB.  At every target, a separate covariance-domain
% communication carrier is also propagated and competes with the best-CRB
% state only after a complete AO solve at the new target.
dualState.enable = true;
dualState.windowByL = cell(numel(L_list),1);
for iLayer=1:numel(L_list)
    dualState.windowByL{iLayer}=[-Inf,Inf];
end
dualState.etaList = [0.50,0.25]; % linear-SINR fractions toward next target
dualState.maxOuter = 40;
dualState.maxOuterByL = [0,40,60];
dualState.maxPhaseInner = 5;
dualState.maxPhaseInnerByL = [0,5,6];
dualState.outerMinIterations = 5;
dualState.outerStableIterations = 3;
dualState.tolOuter = 1e-4;
dualState.minRelEig = 1e-11; % compatibility only; protected path uses relative drop
dualState.hardMinRelEig = 1e-11;
dualState.maxEigDropFactor = 10;
dualState.maxCRBStepDB = 3.0;
dualState.minStepDBByL = [Inf,0.03125,0.015625];
dualState.maxProtectedAttempts = 48;
dualState.phaseTrustRadius0ByL = [NaN,0.06,0.04];
dualState.phaseTrustRadiusMaxByL = [NaN,0.12,0.08];
dualState.reachRelativeTolerance = 1e-10;
dualState.maxCRBLossDB = Inf; % transition state only needs feasibility/stability
dualState.verbose = true;

% ----------------------- Persisted transition diagnostics ----------------
% Diagnostics are written to a separate MAT file and do not alter the
% formal result MAT file.  A fixed-schema table is used deliberately: no
% empty structure is assigned to another structure or structure array.
diagnosticOptions.enable = true;
diagnosticOptions.saveAfterEveryAttempt = true;
diagnosticOptions.printBridgeAttempts = true;
diagnosticOptions.writeCompactPointTrace = true;
diagnosticOptions.maxAbnormalDetailsPerLayer = 12;

refine.enable = true;

% Manually monitored jump regions inferred from the current FP-SCA result.
% Every overlapping original 1-dB interval is divided into four equal parts
% in the linear-SINR domain. Hence the monitored density is approximately
% 0.25 dB without refining the whole curve.
refine.forceWindowByL = cell(numel(L_list),1);
refine.forceWindowByL{1} = [];
refine.forceWindowByL{2} = [13,24];
refine.forceWindowByL{3} = [10.5,25];
if isRRIS
    refine.forceWindowByL={[],[],[]};
end
refine.forceSubdivisions = 4;
refine.forceWindowDepth = ceil(log2(refine.forceSubdivisions));

% Forced monitoring starts at depth 2. Soft transitions may reach depth 3,
% while hard transitions or failures may reach depth 4.
refine.softMaxDepth = 5;
refine.hardMaxDepth = 6;
refine.maxDepth = refine.hardMaxDepth; % compatibility alias
refine.maxInsertedPerGap = 12;
refine.maxTrialsPerGap = 28;
refine.minRelativeGammaGap = 5e-4;
refine.softCRBJumpDB = 1.0;
refine.softEigDropFactor = 5;
refine.hardCRBJumpDB = 3.0;
refine.hardEigDropFactor = 20;
refine.phasePerturbStd = 0.05;
refine.seed = params0.seed.init + 180000;

% Compatibility aliases retained in saved MAT files.
abnormalCRBJumpDB = refine.hardCRBJumpDB;
abnormalEigDropDecades = log10(refine.hardEigDropFactor);

% Reproducible multi-start initialization, strongest for L=3.
rangeInitOpts.profile = 'multistart';
rangeInitOpts.seed = params0.seed.init + 90000;
rangeInitOpts.maxCandidates = 3; % overwritten by L: 3/5/8
rangeInitOpts.randomizeW = false;
rangeInitCountByL = [3,5,8];

anchorInitOpts.profile = 'multistart';
anchorInitOpts.seed = params0.seed.init + 120000;
anchorInitOpts.maxCandidates = 8;
anchorCandidateCapByL = [2,4,8];
anchorInitOpts.randomizeW = false;

recoveryInitOpts.profile = 'compact';
recoveryInitOpts.seed = params0.seed.init + 150000;
recoveryInitOpts.maxCandidates = 2;
recoveryInitOpts.randomizeW = false;

% ---------------- Pure-communication boundary / high-SINR rescue ---------
% If an independently verified pure-communication limit is already known,
% put it here in dB. Example: commBoundary.referenceDBByL(3) = 29;
% NaN means: estimate it automatically with a reproducible communication-
% only multi-start fixed-phase SDR search.
commBoundary.enable = true;
commBoundary.referenceDBByL = [NaN,23.1,24.3];
if isRRIS
    commBoundary.referenceDBByL=[NaN,NaN,NaN];
end
commBoundary.safetyMarginDB = 0.05;
commBoundary.seed = params0.seed.init + 260000;
commBoundary.numRandomByL = [0,6,10];
commBoundary.maxCandidatesByL = [1,10,14];
commBoundary.gammaSeed = 1;
commBoundary.maxExpand = 35;
commBoundary.maxBisection = 32;
commBoundary.upperCapDB = 45;
commBoundary.relativeTolerance = 2e-4;
commBoundary.cvxQuiet = true;
commBoundary.cvxPrecisionMode = 'default';
commBoundary.nearBoundaryWidthDB = 5.0;
commBoundary.midBoundaryWidthDB = 2.0;
commBoundary.veryNearBoundaryWidthDB = 0.5;
commBoundary.nearBoundaryStepDB = 0.25;
commBoundary.midBoundaryStepDB = 0.10;
commBoundary.veryNearBoundaryStepDB = 0.025;
commBoundary.rescueMaxOuterByL = [0,300,360];
commBoundary.rescueMaxPhaseInnerByL = [0,8,10];
commBoundary.rescueNumPerturbByL = [0,2,3];
commBoundary.rescuePerturbStdByL = [0,0.03,0.02];
commBoundary.fixedFallbackRegRel = 1e-7;

% --------------------------- AO parameters -------------------------------
crbAlg = params0.alg;

% Balanced test profile. The outer caps are lower than the old fixed-budget
% AL run because every FP-SCA phase iteration solves a small convex model.
% Early stopping prevents stable points from consuming the full cap.
crbAlg.phaseMethod = 'FP-SCA';
crbAlg.maxOuter = 120;
crbAlg.maxPhaseInner = 5;
crbAlg.useOuterStop = true;
crbAlg.tolOuter = 1e-5;
crbAlg.outerMinIterations = 15;
crbAlg.outerStableIterations = 6;
crbAlg.tolPhase = 0;
crbAlg.phaseGradTol = 1e-7;
crbAlg.debugBlock = false;
crbAlg.verbose = false;

% The physical CRB is never loaded or clipped. The warning threshold is
% diagnostic; the hard threshold rejects numerically unreliable candidates.
crbAlg.crbWarnMinRelEig = 1e-8;
crbAlg.crbHardMinRelEig = 1e-12;
crbAlg.crbMinRcond = crbAlg.crbHardMinRelEig; % compatibility alias
crbAlg.phaseGuardMinRelEig = 1e-8;
crbAlg.phaseMaxRelativeEigDrop = 5;
crbAlg.useLogCRBPhaseObjective = true;
% The phase direction, every acceptance decision, and the saved curve use
% the same unregularized physical CRB.  The fixed-phase SDR retains its
% adaptive numerical loading because a no-loading trial is infeasible at all
% three diagnosed L=2/L=3 branch transitions.  Every returned SDR candidate
% is still checked with the exact raw CRB and covariance-domain SINR.
crbAlg.useInternalCRBRegularization = false;
crbAlg.crbInternalRegRel = 0; % compatibility alias for phase gradient
crbAlg.crbInternalRegRelA = 0;
crbAlg.crbInternalRegRelB = 0;
crbAlg.useRobustSDR = true;
crbAlg.robustSDRAdaptive = true;
crbAlg.robustSDRSchurScale = 'N';
crbAlg.robustSDRFixedFallbackRegRel = 1e-7;
crbAlg.robustSDRFinalBestRetry = false;
crbAlg.cvxPrecisionMode = 'high';
crbAlg.cvxSchurScale = 1;
if isRRIS
    % Multi-hop attenuation is removed analytically by whitening/scaling.
    % Use tighter solver accuracy for the residual conditioning problem;
    % acceptance still uses the raw, unregularized CRB and original SINR.
    crbAlg.robustSDRPrecisionMode='high';
    crbAlg.robustSDRFinalBestRetry=true;
    crbAlg.robustSDRFinalPrecisionMode='best';
end
% Common loading used only for the diagnostic regularized-CRB figure.
% The MAT file stores only exact unregularized CRB values.
crbAlg.crbPlotRegRel = 0;
crbAlg.sensingEigFloorRel = 1e-8;

% Feasible monotone proximal-SCA phase parameters. The exact CRB and exact
% scaled SINR inequalities are checked after every trial. A direction is
% solved at most twice per phase step, while inexpensive line-search trials
% are used first. Each AO outer iteration now permits 5 phase steps in the
% ordinary region and 6 in the manually refined region.
crbAlg.phaseConstraintTol = 2e-5;
crbAlg.phaseExactConstraintTol = 2e-5;
crbAlg.phaseSubproblemConstraintTol = 2e-5;
crbAlg.phaseTrustRadius0 = 0.08;
crbAlg.phaseTrustRadiusMin = 1e-5;
crbAlg.phaseTrustRadiusMax = 0.20;
crbAlg.phaseTrustShrink = 0.5;
crbAlg.phaseTrustExpand = 1.10;
crbAlg.phaseCurvatureF0 = 1;
crbAlg.phaseCurvatureC0 = 1;
crbAlg.phaseCurvatureGrowth = 2;
crbAlg.phaseCurvatureDecrease = 1.10;
crbAlg.phaseCurvatureMin = 1e-4;
crbAlg.phaseCurvatureMax = 1e8;
crbAlg.phaseMaxModelResolve = 2;
crbAlg.phaseMaxLineSearch = 8;
crbAlg.phaseLineSearchBeta = 0.5;
crbAlg.phaseMinStepScale = 1/128;
crbAlg.phaseStepTol = 1e-6;
crbAlg.phaseModelTol = 1e-7;
crbAlg.phaseObjectiveTol = 1e-9;
crbAlg.phaseSufficientDecrease = 1e-4;
% The local quadratic model generates a direction but is not required to be
% a certified upper bound.  Exact SINR feasibility and exact physical-CRB
% monotonicity below remain the hard acceptance tests.
crbAlg.phaseRequireModelUpperBound = false;

% The CRB-invariant outermost phase Phi_L is optimized only for true
% covariance-domain SINR reserve. The inner CRB block then retains a small
% fraction of that reserve and updates only Phi_1,...,Phi_{L-1}. For L=1
% there are no inner sensing phases, so the CRB phase block is skipped.
crbAlg.outerMarginEnable = true;
crbAlg.outerMarginMaxIter = 12;
crbAlg.outerMarginMaxLineSearch = 10;
crbAlg.outerMarginAlpha0 = 0.20;
crbAlg.outerMarginAlphaMin = 1e-5;
crbAlg.outerMarginLineSearchBeta = 0.5;
crbAlg.outerMarginSmoothTau = 100;
crbAlg.outerMarginGradTol = 1e-8;
crbAlg.outerMarginImproveTol = 1e-7;
crbAlg.outerMarginCRBRelTol = 1e-9;
crbAlg.innerReserveFraction = 0.10;
crbAlg.innerReserveFloor = 1e-5;

crbAlg.cvxQuiet = true;
crbAlg.progressEvery = 25;
crbAlg.cvxPowerTol = 1e-8;
crbAlg.sdrPostSINRRelTol = 1e-5;
crbAlg.sdrPostSINRAbsTol = 1e-10;
crbAlg.recoveryEigTol = 1e-12; % compatibility alias
crbAlg.recoveryDenomRelTol = 1e-12;
crbAlg.recoveryPSDSoftTol = 1e-10;
crbAlg.recoveryPSDHardTol = 1e-7;
crbAlg.recoveryCovarianceGapTol = 1e-7;
crbAlg.recoveryDesiredGapTol = 1e-8;
crbAlg.recoverySelfLeakTol = 1e-8;
crbAlg.blockAcceptRelTol = 1e-9;
crbAlg.maxSDRRetry = 2;
crbAlg.maxConsecutiveSDRFailures = 3;

% Persistently propagated communication carrier.  The baseline reserve is
% 10%% of the CURRENT linear SINR, matching the reference projects, and the
% target grows with the actual next-point step up to a 50%% ceiling.
crbAlg.carrierCurrentTargetMarginRel = 0.10;
crbAlg.carrierMaxCurrentTargetMarginRel = 0.50;
crbAlg.carrierStepReserveFactor = 1.25;
crbAlg.carrierReserveSatisfactionRatio = 0.85;
crbAlg.carrierFallbackDeliveryMarginRel = 0.01;
crbAlg.carrierSDRSINRGuardRel = 5e-4;
crbAlg.carrierSINRRelativeTol = 2e-5;
crbAlg.carrierFormalMarginFloor = 0;
crbAlg.carrierNumPhaseSeeds = 1;
crbAlg.carrierPhaseIMaxIter = 32;
crbAlg.carrierAllPhaseMaxIter = 50;
crbAlg.carrierOuterMaxIter = 30;
crbAlg.carrierGammaMinStepDB = 1e-3;
crbAlg.carrierCovarianceFloorRelList = [1e-3,1e-5,0];
crbAlg.carrierRequireCRBHealth = true;
crbAlg.carrierMaxCRBEigDropFactor = 10;
crbAlg.carrierMaxCRBLossDB = 6;
crbAlg.carrierRelaxedMaxCRBLossDB = 15;
crbAlg.carrierBlendWeights = logspace(-6,0,25);

% Only the manually refined SINR windows use the larger budget. All range
% searches, anchors, and formal points outside those windows keep crbAlg.
criticalAlg = crbAlg;
criticalAlg.maxOuter = 240;
criticalAlg.maxPhaseInner = 8;
criticalAlg.outerMinIterations = 20;
criticalAlg.outerStableIterations = 8;
criticalAlg.crbInternalRegRel = 0;
criticalAlg.crbInternalRegRelA = 0;
criticalAlg.crbInternalRegRelB = 0;
criticalAlg.useRobustSDR = true;
criticalAlg.robustSDRAdaptive = true;
criticalAlg.sensingEigFloorRel = 1e-8;
criticalAlg.phaseGuardMinRelEig = 1e-8;
criticalAlg.phaseMaxRelativeEigDrop = 5;
criticalAlg.outerMarginMaxIter = 20;

resultsDir = fullfile(packageDir,'results');
if exist(resultsDir,'dir') ~= 7
    mkdir(resultsDir);
end
fileTag=upper(scenarioCRB);
checkpointFile = fullfile(resultsDir,['CRB_SINR_' fileTag '_checkpoint.mat']);
finalFile = fullfile(resultsDir,['CRB_SINR_' fileTag '_results.mat']);
diagnosticsFile = fullfile(resultsDir,[fileTag '_CRB_SINR_diagnostics.mat']);
diagnosticOptions.traceTextFile = fullfile( ...
    resultsDir,[fileTag '_CRB_SINR_point_diagnostics.txt']);

% A completed layer is expensive and is already self-contained in the
% checkpoint.  Only rRIS enables automatic resume: SIM retains its validated
% fresh-run behavior.  The checkpoint is accepted only when its scenario and
% physical parameters exactly match the current run.
resumeRRIS = false;
resumeStartIndex = 1;
simCRBResults = [];
gammaGridByL = cell(numel(L_list),1);
transitionDiagnostics = make_empty_transition_diagnostic_table();
if isRRIS && exist(checkpointFile,'file') == 2
    savedRun = load(checkpointFile,'simCRBResults','params0', ...
        'gammaGridByL','simCommon','scenarioCRB');
    validSavedRun = isfield(savedRun,'simCRBResults') && ...
        ~isempty(savedRun.simCRBResults) && isfield(savedRun,'params0') && ...
        isequaln(savedRun.params0,params0) && ...
        isfield(savedRun,'scenarioCRB') && strcmpi(savedRun.scenarioCRB,'rRIS');
    if validSavedRun
        savedLayers = [savedRun.simCRBResults.L];
        validSavedRun = isequal(savedLayers,1:numel(savedLayers)) && ...
            numel(savedLayers) < numel(L_list);
    end
    if validSavedRun
        resumeRRIS = true;
        simCRBResults = savedRun.simCRBResults;
        gammaGridByL = savedRun.gammaGridByL;
        simCommon = savedRun.simCommon;
        resumeStartIndex = numel(simCRBResults)+1;
        if exist(diagnosticsFile,'file') == 2
            savedDiag = load(diagnosticsFile,'transitionDiagnostics');
            if isfield(savedDiag,'transitionDiagnostics')
                transitionDiagnostics = savedDiag.transitionDiagnostics;
            end
        end
        diagnosticOptions.traceTextFile = fullfile(resultsDir, ...
            [fileTag '_CRB_SINR_point_diagnostics_resume.txt']);
        initialize_compact_transition_trace(diagnosticOptions,scriptVersion);
        fprintf(['[%s] Resuming from verified checkpoint: completed ', ...
            'L=1:%d; next layer L=%d.\n'],scenarioCRB, ...
            numel(simCRBResults),L_list(resumeStartIndex));
    end
end
if ~resumeRRIS
    initialize_compact_transition_trace(diagnosticOptions,scriptVersion);
    save(diagnosticsFile,'transitionDiagnostics','scriptVersion','dualState', ...
        'refine','diagnosticOptions','-v7.3');
end

fprintf('\n==================== CRB-SINR: %s ====================\n',scenarioCRB);
if resumeRRIS
    % Reuse the exact common channel realization stored with completed layers.
elseif isRRIS
    simCommon=rRIS_build_common_bank(params0,max(L_list));
else
    simCommon=SIM_build_common_links(params0);
end

Lref_SIM = 1;
if isRRIS
    chNoiseRefSIM=rRIS_build_channels(params0,Lref_SIM,simCommon);
else
    chNoiseRefSIM=SIM_build_channels(params0,Lref_SIM,modeSIM,simCommon);
end
if isRRIS, noiseSeed=params0.seed.init+30000; else, noiseSeed=params0.seed.init+20000; end
[sigma_c2_SIM,sigma_s2_est,noiseInfoSIM] = estimate_reference_noise_average( ...
    params0,chNoiseRefSIM,noiseSeed,numNoiseAvg);
if isRRIS, sigma_s2_SIM=sigma_s2_est; else, sigma_s2_SIM=1e-15; end

fprintf('[%s] Fixed noise: sigma_c2=%.3e, sigma_s2=%.3e, T=%d\n',scenarioCRB, ...
    sigma_c2_SIM,sigma_s2_SIM,params0.T);
fprintf(['[%s] CRB thresholds: warning %.1e, hard %.1e, ', ...
    'phase guard %.1e\n'],scenarioCRB,crbAlg.crbWarnMinRelEig, ...
    crbAlg.crbHardMinRelEig,crbAlg.phaseGuardMinRelEig);

previousLayerMaxSol = [];
previousLayerGammaMax = NaN;
previousLayerAnchorSol = [];
previousLayerTrade = [];
if resumeRRIS
    lastSaved = simCRBResults(end);
    previousLayerGammaMax = lastSaved.gammaMaxFound;
    previousLayerTrade = lastSaved.trade;
    previousLayerAnchorSol = trade_record_to_solution(lastSaved.trade(1));
    previousLayerMaxSol = struct('theta',lastSaved.trade(end).theta);
    if isfield(lastSaved,'gammaSearchInfo') && ...
            isfield(lastSaved.gammaSearchInfo,'bestTheta') && ...
            ~isempty(lastSaved.gammaSearchInfo.bestTheta)
        previousLayerMaxSol.theta = lastSaved.gammaSearchInfo.bestTheta;
    end
end

for idx = resumeStartIndex:numel(L_list)
    L = L_list(idx);
    fprintf('\n====================== %s L=%d ======================\n',scenarioCRB,L);
    forcedWindowDB = refine.forceWindowByL{L};
    if isempty(forcedWindowDB)
        fprintf('[%s L=%d] no manually refined SINR window.\n',scenarioCRB,L);
    else
        fprintf(['[%s L=%d] manually refined SINR window: ', ...
            '[%.1f, %.1f] dB; AO/phase budget %d/%d.\n'], ...
            scenarioCRB,L,forcedWindowDB(1),forcedWindowDB(2), ...
            criticalAlg.maxOuter,criticalAlg.maxPhaseInner);
    end

    if isRRIS
        ch0=rRIS_build_channels(params0,L,simCommon);
    else
        ch0=SIM_build_channels(params0,L,modeSIM,simCommon);
    end
    % All layer solvers must see the same reference-calibrated noise.  Most
    % AO routines receive it through initBase, while the independent pure-
    % communication boundary estimator reads it directly from the channel.
    ch0.sigma_c2 = sigma_c2_SIM;
    ch0.sigma_s2 = sigma_s2_SIM;
    initBase = initialize_solver_state(params0,ch0,params0.seed.init);
    initBase.sigma_c2 = sigma_c2_SIM;
    initBase.sigma_s2 = sigma_s2_SIM;

    % ---------------------------------------------------------------------
    % Range search: choose initialization by the largest verified SINR.
    % ---------------------------------------------------------------------
    initOptsL = rangeInitOpts;
    initOptsL.seed = rangeInitOpts.seed + 1000*L;
    initOptsL.maxCandidates = rangeInitCountByL(L);
    [rangeInitList,rangeInitLabels] = SIM_build_cross_layer_inits( ...
        params0,ch0,initBase,previousLayerMaxSol,initOptsL);

    if L == 1
        gammaRequestedStart = gammaTradeStartSIM(1);
        searchOptsL = searchOptsMultiplicative;
        searchOptsL.recoveryMinGamma = gammaTradeStartSIM(1);
    elseif isRRIS
        % Each additional rRIS hop lowers the communication ceiling. Start
        % from this layer's conservative target and expand upward instead of
        % inheriting the preceding layer's much larger maximum SINR.
        gammaRequestedStart = gammaTradeStartSIM(L);
        searchOptsL = searchOptsMultiplicative;
        searchOptsL.maxOuterGrowth = searchOptsAdditive.maxOuterGrowth;
        searchOptsL.maxOuterConfirm = searchOptsAdditive.maxOuterConfirm;
        searchOptsL.maxOuterRefine = searchOptsAdditive.maxOuterRefine;
        searchOptsL.maxExpand = searchOptsAdditive.maxExpand;
        searchOptsL.enableRecovery = true;
        searchOptsL.recoveryFactor = 0.5;
        searchOptsL.maxRecovery = 30;
        searchOptsL.recoveryMinGamma = 1e-8;
    else
        gammaRequestedStart = previousLayerGammaMax;
        searchOptsL = searchOptsAdditive;
        searchOptsL.recoveryMinGamma = gammaTradeStartSIM(L);
        searchOptsL.additiveStep = max(additiveStepAbsMin, ...
            additiveStepRel*previousLayerGammaMax);
        searchOptsL.additiveStepMin = max(additiveStepMinAbs, ...
            additiveStepMinRel*previousLayerGammaMax);
    end

    fprintf(['[%s L=%d] range start %.3f dB with %d initialization(s)\n'], ...
        scenarioCRB,L,10*log10(gammaRequestedStart),numel(rangeInitList));

    rangeTrials = [];
    bestRangeFound = false;
    gammaMaxFound = NaN;
    solAtGammaMax = [];
    searchInfo = [];
    selectedRangeInit = initBase;
    selectedRangeInitLabel = 'none';

    for ii = 1:numel(rangeInitList)
        fprintf('  [range %02d/%02d] %s\n', ...
            ii,numel(rangeInitList),rangeInitLabels{ii});

        searchOptsTrial = searchOptsL;
        searchOptsTrial.enableFullBudgetBoundary = false;
        [gammaTmp,solTmp,infoTmp] = SIM_search_max_common_sinr_AO( ...
            params0,ch0,crbAlg,rangeInitList{ii}, ...
            gammaRequestedStart,searchOptsTrial);

        rangeTrials(ii).label = rangeInitLabels{ii};
        rangeTrials(ii).success = isfield(infoTmp,'success') && ...
            logical(infoTmp.success) && isfinite(gammaTmp);
        rangeTrials(ii).gammaMaxFound = gammaTmp;
        rangeTrials(ii).gammaMaxFoundDB = ...
            10*log10(max(gammaTmp,realmin));
        rangeTrials(ii).searchInfo = infoTmp;
        rangeTrials(ii).bestCRB = candidate_crb(solTmp);

        if ~rangeTrials(ii).success
            continue;
        end

        replaceBest = ~bestRangeFound || gammaTmp > gammaMaxFound*(1+1e-8);
        if bestRangeFound && abs(gammaTmp-gammaMaxFound) <= ...
                1e-8*max(1,gammaMaxFound)
            replaceBest = candidate_crb(solTmp) < candidate_crb(solAtGammaMax);
        end

        if replaceBest
            bestRangeFound = true;
            gammaMaxFound = gammaTmp;
            solAtGammaMax = solTmp;
            searchInfo = infoTmp;
            selectedRangeInit = rangeInitList{ii};
            selectedRangeInitLabel = rangeInitLabels{ii};
        end
    end

    if ~bestRangeFound
        error('[%s L=%d] No initialization produced a feasible SINR range.',scenarioCRB,L);
    end

    fprintf('[%s L=%d] boundary confirmation from %s\n', ...
        scenarioCRB,L,selectedRangeInitLabel);
    [gammaConfirmed,solConfirmed,infoConfirmed] = ...
        SIM_search_max_common_sinr_AO(params0,ch0,crbAlg, ...
        selectedRangeInit,gammaRequestedStart,searchOptsL);

    if isfield(infoConfirmed,'success') && infoConfirmed.success && ...
            isfinite(gammaConfirmed) && ...
            gammaConfirmed >= gammaMaxFound*(1-1e-8)
        gammaMaxFound = gammaConfirmed;
        solAtGammaMax = solConfirmed;
        searchInfo = infoConfirmed;
    else
        warning('SIM_CRB_SINR:BoundaryConfirmationFailed', ...
            ['%s L=%d full-budget boundary confirmation failed or ', ...
             'returned a lower range. The best verified short-budget ', ...
             'range is retained.'],scenarioCRB,L);
    end

    fprintf('[%s L=%d] largest CRB-AO verified Gamma %.3f dB\n', ...
        scenarioCRB,L,10*log10(gammaMaxFound));

    % The CRB-AO range is not the communication-side physical limit. For
    % L=2/L=3 independently estimate (or manually supply) the pure-
    % communication common-SINR ceiling and let the formal CRB curve continue
    % toward that ceiling. The exact endpoint is backed off by a tiny safety
    % margin because finite raw CRB may diverge at the communication optimum.
    if L >= 2 && commBoundary.enable
        [gammaCommMax,commBoundaryInfo] = ...
            SIM_estimate_pure_comm_sinr_boundary( ...
            params0,ch0,rangeInitList,commBoundary);
        if commBoundaryInfo.success
            gammaFormalMax = max(gammaMaxFound, ...
                gammaCommMax*10^(-commBoundary.safetyMarginDB/10));
            fprintf(['[%s L=%d] pure-communication ceiling %.3f dB; ', ...
                'formal CRB sweep endpoint %.3f dB.\n'],scenarioCRB,L, ...
                commBoundaryInfo.gammaCommMaxDB,10*log10(gammaFormalMax));
        else
            gammaFormalMax = gammaMaxFound;
            warning('SIM_CRB_SINR:PureCommBoundaryFailed', ...
                ['%s L=%d pure-communication boundary estimate failed; ', ...
                 'falling back to the CRB-AO verified endpoint.'],scenarioCRB,L);
        end
    else
        gammaFormalMax = gammaMaxFound;
        commBoundaryInfo = struct( ...
            'success',false,'gammaCommMax',NaN,'gammaCommMaxDB',NaN, ...
            'bestTheta',zeros(0,0),'bestLabel','none', ...
            'status','L1 uses original CRB-AO endpoint');
    end

    if L > 1 && gammaMaxFound < previousLayerGammaMax
        if isRRIS
            fprintf(['[rRIS L=%d] verified ceiling %.3f dB is below ', ...
                'L=%d ceiling %.3f dB, as expected for an additional hop.\n'], ...
                L,10*log10(gammaMaxFound),L-1,10*log10(previousLayerGammaMax));
        else
            warning('SIM_CRB_SINR:NonMonotoneRange', ...
                ['SIM L=%d verified range remains below L=%d: %.3f dB ', ...
                 'versus %.3f dB.'],L,L-1,10*log10(gammaMaxFound), ...
                 10*log10(previousLayerGammaMax));
        end
    end

    gammaGridStart = min(gammaTradeStartSIM(L),gammaFormalMax);
    if L >= 2 && gammaFormalMax > gammaMaxFound*(1+1e-10)
        gammaGridCoarse = SIM_build_boundary_aware_gamma_grid( ...
            gammaGridStart,gammaFormalMax,formalStepDB,commBoundary);
    else
        gammaGridCoarse = build_uniform_db_grid( ...
            gammaGridStart,gammaFormalMax,formalStepDB);
    end
    gammaGrid = gammaGridCoarse;
    gammaGridByL{idx} = gammaGrid;

    fprintf(['[%s L=%d] coarse formal grid %.3f to %.3f dB, ', ...
        '%d point(s), step <= %.2f dB\n'],scenarioCRB,L, ...
        10*log10(gammaGridCoarse(1)), ...
        10*log10(gammaGridCoarse(end)), ...
        numel(gammaGridCoarse),formalStepDB);

    % The retained grid is dynamic because successful local-refinement
    % points are inserted as complete formal records.
    trade = empty_trade_record();

    % ---------------------------------------------------------------------
    % First formal point: CRB-oriented multi-initialization anchor.
    % ---------------------------------------------------------------------
    anchorOptsL = anchorInitOpts;
    anchorOptsL.seed = anchorInitOpts.seed + 1000*L;
    [compactAnchorList,compactAnchorLabels] = ...
        SIM_build_cross_layer_inits( ...
        params0,ch0,initBase,previousLayerAnchorSol,anchorOptsL);

    % Anchor candidates: keep the range-selected start and add
    % reproducible independent alternatives. L=3 may retain up to eight
    % candidates because its low-SINR branch showed initialization sensitivity.
    anchorInitList = {};
    anchorLabels = {};
    [anchorInitList,anchorLabels] = append_unique_init( ...
        anchorInitList,anchorLabels,selectedRangeInit, ...
        ['selected-range-',selectedRangeInitLabel]);
    preferredAnchorOrder = find( ...
        contains(compactAnchorLabels,'structured') | ...
        contains(compactAnchorLabels,'zero'));
    remainingAnchorOrder = setdiff(1:numel(compactAnchorList), ...
        preferredAnchorOrder,'stable');
    preferredAnchorOrder = [preferredAnchorOrder(:); ...
        remainingAnchorOrder(:)];
    for ia = reshape(preferredAnchorOrder,1,[])
        [anchorInitList,anchorLabels] = append_unique_init( ...
            anchorInitList,anchorLabels,compactAnchorList{ia}, ...
            compactAnchorLabels{ia});
        if numel(anchorInitList) >= anchorCandidateCapByL(L)
            break;
        end
    end
    if numel(anchorInitList) < 2
        [anchorInitList,anchorLabels] = append_unique_init( ...
            anchorInitList,anchorLabels,initBase,'current-layer-base');
    end

    gammaValue = gammaGridCoarse(1);
    gammaTarget = gammaValue*ones(params0.K,1);
    fprintf('\n  [%s L=%d anchor] target %.3f dB, %d initialization(s)\n', ...
        scenarioCRB,L,10*log10(gammaValue),numel(anchorInitList));

    anchorCandidates = [];
    for ii = 1:numel(anchorInitList)
        [solTry,histTry] = SIM_run_CRB_SINR_AO_solver( ...
            params0,ch0,crbAlg,anchorInitList{ii},gammaTarget);
        anchorCandidates = append_solution_candidate( ...
            anchorCandidates,anchorLabels{ii},solTry,histTry,gammaValue);
        print_candidate_summary(anchorLabels{ii},solTry,params0);
    end

    [anchorBest,anchorFound] = select_best_complete_candidate(anchorCandidates);
    if ~anchorFound
        error('[%s L=%d] No valid CRB-oriented formal anchor was found.',scenarioCRB,L);
    end

    anchorSummary = candidates_to_summary(anchorCandidates);
    trade(1) = make_trade_record(anchorBest.sol,anchorBest.hist, ...
        gammaValue,gammaTarget,gammaGridStart,gammaMaxFound,params0, ...
        anchorBest.label,anchorSummary,false,'');
    trade(1).useCriticalBudget = false;
    trade(1).maxOuterUsed = crbAlg.maxOuter;
    trade(1).maxPhaseInnerUsed = crbAlg.maxPhaseInner;
    warmInit = solution_to_init(initBase,anchorBest.sol, ...
        sigma_c2_SIM,sigma_s2_SIM);
    commCarrierWarmInit=warmInit;
    commCarrierHistory={commCarrierWarmInit};
    commCarrierDebugHistory={struct('status','anchor-bestCRB', ...
        'nextMargin',NaN,'targetCurrentMargin',NaN,'targetNextMargin',NaN, ...
        'crbDeltaDB',0,'quality','protected','selected',false)};

    fprintf(['  [anchor selected] %s: CRB %.3f dB, minSINR %.3f dB, ', ...
        'relEig(A/B)=%.2e/%.2e\n'],anchorBest.label, ...
        trade(1).CRBPhysicalDB,trade(1).minSINRdB, ...
        trade(1).relMinEigA,trade(1).relMinEigB);

    anchorTransitionDiag = make_default_transition_diag( ...
        trade(1).targetGammaScalar,trade(1).targetGammaScalar, ...
        'anchor-report');
    anchorTransitionDiag.carrierInfo=struct( ...
        'primaryPath',summarize_ao_path('anchor-bestCRB', ...
        anchorBest.sol,anchorBest.hist,trade(1).targetGammaScalar, ...
        params0.T,0,crbAlg.maxOuter), ...
        'additionalCandidates',summarize_additional_candidate_paths( ...
        anchorCandidates,trade(1).targetGammaScalar,params0.T), ...
        'quality','anchor','buildElapsedSeconds',0);
    transitionDiagnostics = append_transition_diagnostic_row( ...
        transitionDiagnostics,L,1,trade(1),[], ...
        make_refinement_task(trade(1).targetGammaScalar,0,false, ...
        trade(1).targetGammaScalar,trade(1).targetGammaScalar,false), ...
        false,anchorTransitionDiag,'none',empty_transition(), ...
        false,NaN,false,trade(1),'anchor',params0.T);
    save_transition_diagnostics(diagnosticsFile,transitionDiagnostics, ...
        scriptVersion,dualState,refine,diagnosticOptions);

    gammaGrid = [trade.targetGammaScalar].';
    gammaGridByL{idx} = gammaGrid;
    resultEntry = make_SIM_result_entry( ...
        L,ch0.mode,ch0,initBase,gammaGridStart,gammaRequestedStart, ...
        gammaMaxFound,gammaGrid,formalStepDB,searchInfo,rangeTrials, ...
        selectedRangeInitLabel,trade,noiseInfoSIM,sigma_c2_SIM, ...
        sigma_s2_SIM,params0.T);
    resultEntry.gammaGridCoarse = gammaGridCoarse;
    resultEntry.refineOptions = refine;
    resultEntry.gammaFormalMax = gammaFormalMax;
    resultEntry.commBoundaryInfo = commBoundaryInfo;
    resultEntry.commCarrier=commCarrierHistory;
    resultEntry.commCarrierDebug=commCarrierDebugHistory;
    simCRBResults = assign_result_entry( ...
        simCRBResults,idx,resultEntry,'simCRBResults');
    save(checkpointFile,'simCRBResults','params0','crbAlg','criticalAlg', ...
        'searchOptsMultiplicative','searchOptsAdditive', ...
        'gammaSeedSIM','gammaGridByL','formalStepDB','refine','dualState','commBoundary', ...
        'abnormalCRBJumpDB','abnormalEigDropDecades', ...
        'backwardPolishRelTol','enableBackwardPolish', ...
    'enableMonotoneEnvelope','monotoneEnvelopeRelTol', ...
        'enableBackwardReoptimization','backwardRepairTolDB', ...
        'maxBackwardRepairPass','scriptVersion','simCommon','scenarioCRB','-v7.3');

    % ---------------------------------------------------------------------
    % Remaining points: retain the 1-dB coarse targets. Manually selected
    % SINR windows are divided into four linear-SINR subintervals, while
    % every interval may still receive adaptive soft/hard refinement.
    % ---------------------------------------------------------------------
    stopFormalLayer = false;
    for coarseIndex = 2:numel(gammaGridCoarse)
        coarseHighGamma = gammaGridCoarse(coarseIndex);
        lowForGap = last_valid_trade_record(trade);
        coarseLowGamma = lowForGap.targetGammaScalar;
        forceThisGap = is_in_forced_window( ...
            L,coarseLowGamma,coarseHighGamma,refine);

        if forceThisGap
            numSub = max(2,round(refine.forceSubdivisions));
            pendingTasks = repmat(make_refinement_task( ...
                coarseHighGamma,refine.forceWindowDepth,false, ...
                coarseLowGamma,coarseHighGamma,true),numSub,1);
            for isub = 1:numSub
                frac = isub/numSub;
                gammaSub = coarseLowGamma + ...
                    frac*(coarseHighGamma-coarseLowGamma);
                pendingTasks(isub) = make_refinement_task( ...
                    gammaSub,refine.forceWindowDepth,isub<numSub, ...
                    coarseLowGamma,coarseHighGamma,true);
            end
            gapInserted = numSub-1;
            fprintf(['\n  [%s L=%d coarse %02d/%02d] monitored ', ...
                'linear-SINR interval divided into %d parts over ', ...
                '[%.3f, %.3f] dB.\n'], ...
                scenarioCRB,L,coarseIndex,numel(gammaGridCoarse),numSub, ...
                10*log10(coarseLowGamma),10*log10(coarseHighGamma));
        else
            pendingTasks = make_refinement_task( ...
                coarseHighGamma,0,false,coarseLowGamma, ...
                coarseHighGamma,false);
            gapInserted = 0;
        end
        gapTrials = 0;

        while ~isempty(pendingTasks)
            task = pendingTasks(1);
            pendingTasks(1) = [];
            previousRecord = last_valid_trade_record(trade);
            reportWarmInit = trade_record_to_init(initBase,previousRecord, ...
                sigma_c2_SIM,sigma_s2_SIM);
            gammaValue = task.gamma;
            gammaTarget = gammaValue*ones(params0.K,1);

            if task.isForcedPoint
                algThisPoint = criticalAlg;
            else
                algThisPoint = crbAlg;
            end

            % Dual-state continuation: the retained previous record remains
            % the report solution, while an optional bridge state is built
            % only to initialize this next target. The protected internal state is neither
            % appended to trade nor saved as a Pareto point.
            useDualStateHere = is_in_dual_state_window( ...
                L,previousRecord.targetGammaScalar,gammaValue,dualState);
            transitionDiag = make_default_transition_diag( ...
                previousRecord.targetGammaScalar,gammaValue,'report-seed');
            % Best-CRB and communication-carrier states are independent and
            % both are propagated across every formal SINR point.
            warmInit = reportWarmInit;

            fprintf(['\n  [%s L=%d coarse %02d/%02d, depth %d] ', ...
                'target %.3f dB%s%s, budget %d/%d\n'],scenarioCRB,L,coarseIndex, ...
                numel(gammaGridCoarse),task.depth, ...
                10*log10(gammaValue),inserted_suffix(task.isInsertedPoint), ...
                forced_suffix(task.isForcedPoint), ...
                algThisPoint.maxOuter,algThisPoint.maxPhaseInner);

            [bestCandidate,trialCandidates,severity,recoveryReason, ...
                transition,carrierInfo,nextCarrierWarmInit] = ...
                solve_compact_formal_target( ...
                params0,ch0,algThisPoint,warmInit,previousRecord, ...
                gammaValue,gammaTarget,refine,coarseIndex,useDualStateHere, ...
                commCarrierWarmInit);
            transitionDiag.carrierInfo = carrierInfo;
            if carrierInfo.attempted
                transitionDiag.attempted = true;
                transitionDiag.usedBridge = carrierInfo.selected;
                transitionDiag.handoffSource = carrierInfo.handoffSource;
                transitionDiag.etaUsed = carrierInfo.blendWeight;
                transitionDiag.bridgeTarget = gammaValue;
                transitionDiag.bridgeTargetDB = 10*log10(gammaValue);
                transitionDiag.bridgeCRBPhysicalDB = ...
                    10*log10(max((sigma_s2_SIM/params0.T)* ...
                    carrierInfo.carrierCRB,realmin));
                transitionDiag.minSINR = gammaValue* ...
                    (1+carrierInfo.nextMargin);
                transitionDiag.minSINRdB = 10*log10(max( ...
                    transitionDiag.minSINR,realmin));
                transitionDiag.relMinEigA = carrierInfo.relMinEigA;
                transitionDiag.relMinEigB = carrierInfo.relMinEigB;
                transitionDiag.status = carrierInfo.status;
                transitionDiag.numAttempts = carrierInfo.numAttempts;
                transitionDiag.attemptEtaList = carrierInfo.blendWeight;
                transitionDiag.attemptTargetDB = 10*log10(gammaValue);
                transitionDiag.attemptAccepted = carrierInfo.selected;
                transitionDiag.attemptCRBPhysicalDB = ...
                    transitionDiag.bridgeCRBPhysicalDB;
                transitionDiag.attemptMinSINRdB = transitionDiag.minSINRdB;
                transitionDiag.attemptRelMinEigA = carrierInfo.relMinEigA;
                transitionDiag.attemptRelMinEigB = carrierInfo.relMinEigB;
                transitionDiag.attemptStatus = {carrierInfo.status};
            end
            gapTrials = gapTrials + 1;

            relativeGap = (gammaValue-previousRecord.targetGammaScalar) / ...
                max(previousRecord.targetGammaScalar,realmin);
            needsProtection = isempty(bestCandidate) || ...
                ~strcmp(severity,'none');
            switch severity
                case 'soft'
                    allowedDepth = refine.softMaxDepth;
                case {'hard','failure'}
                    allowedDepth = refine.hardMaxDepth;
                otherwise
                    allowedDepth = 0;
            end
            canRefine = refine.enable && ...
                task.depth < allowedDepth && ...
                gapInserted < refine.maxInsertedPerGap && ...
                gapTrials < refine.maxTrialsPerGap && ...
                relativeGap > refine.minRelativeGammaGap;

            if needsProtection && canRefine
                gammaMid = 0.5*(previousRecord.targetGammaScalar+gammaValue);
                retryTask = task;
                retryTask.depth = task.depth+1;
                midTask = make_refinement_task( ...
                    gammaMid,task.depth+1,true, ...
                    previousRecord.targetGammaScalar,gammaValue, ...
                    task.isForcedPoint);
                pendingTasks = [midTask;retryTask;pendingTasks]; %#ok<AGROW>
                gapInserted = gapInserted + 1;
                fprintf(['    local refinement inserted %.3f dB before ', ...
                    'retesting %.3f dB (%s).\n'], ...
                    10*log10(gammaMid),10*log10(gammaValue), ...
                    recoveryReason);
                transitionDiagnostics = append_transition_diagnostic_row( ...
                    transitionDiagnostics,L,coarseIndex,previousRecord, ...
                    bestCandidate,task,useDualStateHere,transitionDiag, ...
                    severity,transition,true,gammaMid,false,[], ...
                    recoveryReason,params0.T);
                save_transition_diagnostics(diagnosticsFile, ...
                    transitionDiagnostics,scriptVersion,dualState,refine, ...
                    diagnosticOptions);
                continue;
            end

            branchProtectionExhausted = needsProtection && ~canRefine;
            branchSwitchTried = false;
            % A feasible continuation is NEVER replaced merely because its
            % CRB jump is large. Large positive jumps are handled only by
            % target refinement. Random/structured branch recovery is allowed
            % after the refinement budget is exhausted only when continuation
            % itself could not produce a feasible complete solution.
            if isempty(bestCandidate)
                branchSwitchTried = true;
                [branchBest,branchCandidates] = solve_branch_switch_target( ...
                    params0,ch0,algThisPoint,initBase,warmInit, ...
                    previousLayerTrade,gammaValue,gammaTarget, ...
                    recoveryInitOpts,refine,coarseIndex,task.depth);
                trialCandidates = concatenate_solution_candidates( ...
                    trialCandidates,branchCandidates);
                if ~isempty(branchBest)
                    bestCandidate = branchBest;
                end
            end

            % Before terminating L=2/L=3, run a stronger boundary rescue.
            % It uses the last valid report state, the pure-communication
            % phase seed (if available), deterministic small perturbations,
            % and adaptive A/B regularization. Only exact raw-CRB/SINR
            % feasible complete solutions may be returned.
            if isempty(bestCandidate) && L >= 2 && commBoundary.enable
                [rescueCandidate,rescueInfo] = SIM_rescue_high_sinr_target( ...
                    params0,ch0,algThisPoint,initBase,previousRecord, ...
                    gammaValue,gammaTarget,commBoundaryInfo,commBoundary, ...
                    sigma_c2_SIM,sigma_s2_SIM);
                fprintf('    high-SINR rescue: %s\n',rescueInfo.status);
                if rescueInfo.success
                    trialCandidates = concatenate_solution_candidates( ...
                        trialCandidates,rescueCandidate);
                    bestCandidate = rescueCandidate;
                    recoveryReason = [recoveryReason,'; high-SINR-rescue'];
                end
            end
            transitionDiag.carrierInfo.additionalCandidates = ...
                summarize_additional_candidate_paths( ...
                trialCandidates,gammaValue,params0.T);

            if isempty(bestCandidate)
                failSol = point_failure_solution(gammaTarget);
                record = make_trade_record(failSol,[], ...
                    gammaValue,gammaTarget,gammaGridStart,gammaMaxFound, ...
                    params0,'none',candidates_to_summary(trialCandidates), ...
                    true,recoveryReason);
                record.isInsertedPoint = task.isInsertedPoint;
                record.isForcedPoint = task.isForcedPoint;
                record.useCriticalBudget = task.isForcedPoint;
                record.maxOuterUsed = algThisPoint.maxOuter;
                record.maxPhaseInnerUsed = algThisPoint.maxPhaseInner;
                record.refinementDepth = task.depth;
                record.parentLowGamma = task.parentLowGamma;
                record.parentHighGamma = task.parentHighGamma;
                record.branchProtectionExhausted = ...
                    branchProtectionExhausted;
                record.branchSwitchTried = branchSwitchTried;
                trade(end+1,1) = orderfields(record,trade(1)); %#ok<AGROW>
                commCarrierHistory{end+1,1}=[]; %#ok<SAGROW>
                commCarrierDebugHistory{end+1,1}=carrierInfo; %#ok<SAGROW>
                transitionDiagnostics = append_transition_diagnostic_row( ...
                    transitionDiagnostics,L,coarseIndex,previousRecord, ...
                    [],task,useDualStateHere,transitionDiag,'failure', ...
                    empty_transition(),false,NaN,branchSwitchTried,record, ...
                    recoveryReason,params0.T);
                save_transition_diagnostics(diagnosticsFile, ...
                    transitionDiagnostics,scriptVersion,dualState,refine, ...
                    diagnosticOptions);
                fprintf(['    failed after local-refinement/branch-switch ', ...
                    'budget; stopping this L.\n']);
                stopFormalLayer = true;
                break;
            end

            % Re-evaluate transition metadata for the finally selected
            % complete target candidate.
            [severity,recoveryReason,transition] = ...
                classify_solution_transition(bestCandidate.sol, ...
                previousRecord,params0,refine);
            record = make_trade_record(bestCandidate.sol, ...
                bestCandidate.hist,gammaValue,gammaTarget,gammaGridStart, ...
                gammaMaxFound,params0,bestCandidate.label, ...
                candidates_to_summary(trialCandidates),needsProtection, ...
                recoveryReason);
            record.isInsertedPoint = task.isInsertedPoint;
            record.isForcedPoint = task.isForcedPoint;
            record.useCriticalBudget = task.isForcedPoint;
            record.maxOuterUsed = algThisPoint.maxOuter;
            record.maxPhaseInnerUsed = algThisPoint.maxPhaseInner;
            record.refinementDepth = task.depth;
            record.parentLowGamma = task.parentLowGamma;
            record.parentHighGamma = task.parentHighGamma;
            record.transitionSeverity = severity;
            record.deltaCRBDBFromPrevious = transition.deltaCRBDB;
            record.pointDropA = transition.pointDropA;
            record.pointDropB = transition.pointDropB;
            record.branchWarning = ~strcmp(severity,'none');
            record.branchRejected = false;
            record.branchProtectionExhausted = ...
                branchProtectionExhausted;
            record.branchSwitchTried = branchSwitchTried;
            if startsWith(bestCandidate.label,'high-SINR-rescue-')
                record.solutionSource = 'high-SINR-rescue';
            elseif startsWith(bestCandidate.label,'branch-')
                record.solutionSource = 'branch-recovery-after-failure';
            elseif strcmp(bestCandidate.label,'commCarrier')
                record.solutionSource = 'persistent-commCarrier';
            elseif branchProtectionExhausted
                record.solutionSource = ...
                    'continuation-after-refinement-exhausted';
            elseif task.isInsertedPoint
                record.solutionSource = 'local-refinement';
            else
                record.solutionSource = 'direct-solve';
            end

            trade(end+1,1) = orderfields(record,trade(1)); %#ok<AGROW>
            if isempty(nextCarrierWarmInit)
                commCarrierWarmInit=solution_to_init(initBase, ...
                    bestCandidate.sol,sigma_c2_SIM,sigma_s2_SIM);
            else
                commCarrierWarmInit=nextCarrierWarmInit;
            end
            commCarrierHistory{end+1,1}=commCarrierWarmInit; %#ok<SAGROW>
            commCarrierDebugHistory{end+1,1}=carrierInfo; %#ok<SAGROW>

            fprintf(['    selected=%s, source=%s, CRB %.3f dB, ', ...
                'minSINR %.3f dB, relEig(A/B)=%.2e/%.2e, ', ...
                'transition=%s\n'],bestCandidate.label, ...
                record.solutionSource,record.CRBPhysicalDB, ...
                record.minSINRdB,record.relMinEigA,record.relMinEigB, ...
                severity);

            transitionDiagnostics = append_transition_diagnostic_row( ...
                transitionDiagnostics,L,coarseIndex,previousRecord, ...
                bestCandidate,task,useDualStateHere,transitionDiag, ...
                severity,transition,false,NaN,branchSwitchTried,record, ...
                recoveryReason,params0.T);
            save_transition_diagnostics(diagnosticsFile, ...
                transitionDiagnostics,scriptVersion,dualState,refine, ...
                diagnosticOptions);

            gammaGrid = [trade.targetGammaScalar].';
            gammaGridByL{idx} = gammaGrid;
            resultEntry = make_SIM_result_entry( ...
                L,ch0.mode,ch0,initBase,gammaGridStart, ...
                gammaRequestedStart,gammaMaxFound,gammaGrid,formalStepDB, ...
                searchInfo,rangeTrials,selectedRangeInitLabel,trade, ...
                noiseInfoSIM,sigma_c2_SIM,sigma_s2_SIM,params0.T);
            resultEntry.gammaGridCoarse = gammaGridCoarse;
            resultEntry.refineOptions = refine;
    resultEntry.gammaFormalMax = gammaFormalMax;
    resultEntry.commBoundaryInfo = commBoundaryInfo;
            resultEntry.commCarrier=commCarrierHistory;
            resultEntry.commCarrierDebug=commCarrierDebugHistory;
            simCRBResults = assign_result_entry( ...
                simCRBResults,idx,resultEntry,'simCRBResults');
            save(checkpointFile,'simCRBResults','params0','crbAlg','criticalAlg', ...
                'searchOptsMultiplicative','searchOptsAdditive', ...
                'gammaSeedSIM','gammaGridByL','formalStepDB','refine','dualState','commBoundary', ...
                'abnormalCRBJumpDB','abnormalEigDropDecades', ...
                'backwardPolishRelTol','enableBackwardPolish', ...
                'enableMonotoneEnvelope','monotoneEnvelopeRelTol', ...
        'enableBackwardReoptimization','backwardRepairTolDB', ...
        'maxBackwardRepairPass','scriptVersion','simCommon','scenarioCRB','-v7.3');
        end

        if stopFormalLayer
            break;
        end
    end

    gammaGrid = [trade.targetGammaScalar].';
    gammaGridByL{idx} = gammaGrid;

    % Preserve the genuine strict forward continuation path before any
    % cross-target repair.
    tradeForwardRaw = trade;
    [numBackwardViolations,maxBackwardGainDB] = ...
        diagnose_backward_monotonicity( ...
            tradeForwardRaw,monotoneEnvelopeRelTol);
    fprintf(['[%s L=%d] strict forward curve: %d left-bend ', ...
        'violation(s), maximum apparent gain %.3f dB.\n'], ...
        scenarioCRB,L,numBackwardViolations,maxBackwardGainDB);

    % Conditional backward re-optimization. This is NOT a monotone envelope:
    % the higher-target CRB value is never copied to the lower target. The
    % higher-target best W/theta is only used as a valid initialization and
    % incumbent for a genuine lower-target AO solve.
    trade = tradeForwardRaw;
    backwardRepairInfo.enabled = enableBackwardReoptimization;
    backwardRepairInfo.attempts = 0;
    backwardRepairInfo.replacements = 0;
    backwardRepairInfo.maxImprovementDB = 0;
    backwardRepairInfo.rawViolationCount = numBackwardViolations;

    if enableBackwardReoptimization
        for repairPass = 1:maxBackwardRepairPass
            replacementsThisPass = 0;
            for ip = numel(trade)-1:-1:1
                low = trade(ip);
                high = trade(ip+1);
                if ~is_valid_trade_record(low) || ...
                        ~is_valid_trade_record(high)
                    continue;
                end
                if low.CRBPhysicalDB <= ...
                        high.CRBPhysicalDB + backwardRepairTolDB
                    continue;
                end

                if is_in_forced_window( ...
                        L,low.targetGammaScalar,high.targetGammaScalar,refine)
                    algRepair = criticalAlg;
                else
                    algRepair = crbAlg;
                end

                fprintf(['\n  [backward reopt pass %d: %02d <- %02d] ', ...
                    'target %.3f dB, left-bend %.3f dB\n'], ...
                    repairPass,ip,ip+1,low.targetGammaDB, ...
                    low.CRBPhysicalDB-high.CRBPhysicalDB);

                backwardRepairInfo.attempts = ...
                    backwardRepairInfo.attempts + 1;
                [solRepair,histRepair,repairDiag] = ...
                    SIM_backward_reoptimize_CRB_point( ...
                    params0,ch0,algRepair,initBase,high, ...
                    low.targetGammaScalar,low.CRBPhysical, ...
                    backwardRepairTolDB);

                if ~repairDiag.improved
                    fprintf('    no material lower-target improvement retained.\n');
                    continue;
                end

                gammaTargetRepair = ...
                    low.targetGammaScalar*ones(params0.K,1);
                repaired = make_trade_record(solRepair,histRepair, ...
                    low.targetGammaScalar,gammaTargetRepair, ...
                    gammaGridStart,gammaMaxFound,params0, ...
                    sprintf('backward-reopt-from-point-%d',ip+1), ...
                    [],false,'');
                repaired = copy_refinement_metadata(repaired,low);
                repaired.useCriticalBudget = ...
                    algRepair.maxOuter == criticalAlg.maxOuter && ...
                    algRepair.maxPhaseInner == criticalAlg.maxPhaseInner;
                repaired.maxOuterUsed = algRepair.maxOuter;
                repaired.maxPhaseInnerUsed = algRepair.maxPhaseInner;
                repaired.solutionSource = 'backward-reoptimization';
                repaired.backwardReoptimized = true;
                repaired.backwardSourcePoint = ip+1;
                repaired.sourceTargetGamma = high.targetGammaScalar;

                trade(ip) = orderfields(repaired,trade(1));
                replacementsThisPass = replacementsThisPass + 1;
                backwardRepairInfo.replacements = ...
                    backwardRepairInfo.replacements + 1;
                backwardRepairInfo.maxImprovementDB = max( ...
                    backwardRepairInfo.maxImprovementDB, ...
                    repairDiag.improvementDB);
                fprintf(['    retained genuine re-optimization: CRB %.3f -> ', ...
                    '%.3f dB (improvement %.3f dB).\n'], ...
                    low.CRBPhysicalDB,repaired.CRBPhysicalDB, ...
                    repairDiag.improvementDB);
            end
            if replacementsThisPass == 0
                break;
            end
        end
    end

    [finalViolationCount,finalMaxGainDB] = ...
        diagnose_backward_monotonicity(trade,monotoneEnvelopeRelTol);
    backwardRepairInfo.finalViolationCount = finalViolationCount;
    backwardRepairInfo.finalMaxGainDB = finalMaxGainDB;
    fprintf(['[%s L=%d] after backward re-optimization: %d ', ...
        'left-bend violation(s), maximum apparent gain %.3f dB.\n'], ...
        scenarioCRB,L,finalViolationCount,finalMaxGainDB);

    trade = refresh_trade_transition_metadata(trade,refine);

    gammaGrid = [trade.targetGammaScalar].';
    gammaGridByL{idx} = gammaGrid;
    resultEntry = make_SIM_result_entry( ...
        L,ch0.mode,ch0,initBase,gammaGridStart, ...
        gammaRequestedStart,gammaMaxFound,gammaGrid,formalStepDB, ...
        searchInfo,rangeTrials,selectedRangeInitLabel,trade, ...
        noiseInfoSIM,sigma_c2_SIM,sigma_s2_SIM,params0.T);
    resultEntry.gammaGridCoarse = gammaGridCoarse;
    resultEntry.refineOptions = refine;
    resultEntry.gammaFormalMax = gammaFormalMax;
    resultEntry.commBoundaryInfo = commBoundaryInfo;
    resultEntry.commCarrier=commCarrierHistory;
    resultEntry.commCarrierDebug=commCarrierDebugHistory;
    resultEntry.tradeForwardRaw = tradeForwardRaw;
    resultEntry.tradeBackwardReoptimized = [];
    resultEntry.backwardRepairInfo = backwardRepairInfo;
    resultEntry.tradeMonotone = []; % direct monotone envelope removed
    resultEntry.monotoneInfo = struct('enabled',false);
    simCRBResults = assign_result_entry( ...
        simCRBResults,idx,resultEntry,'simCRBResults');

    previousLayerMaxSol = solAtGammaMax;
    previousLayerGammaMax = gammaMaxFound;
    previousLayerTrade = trade;
    previousLayerAnchorSol = trade_record_to_solution(trade(1));

    save(checkpointFile,'simCRBResults','params0','crbAlg','criticalAlg', ...
        'searchOptsMultiplicative','searchOptsAdditive', ...
        'gammaSeedSIM','gammaGridByL','formalStepDB','refine','dualState','commBoundary', ...
        'abnormalCRBJumpDB','abnormalEigDropDecades', ...
        'backwardPolishRelTol','enableBackwardPolish', ...
    'enableMonotoneEnvelope','monotoneEnvelopeRelTol', ...
        'enableBackwardReoptimization','backwardRepairTolDB', ...
        'maxBackwardRepairPass','scriptVersion','simCommon','scenarioCRB','-v7.3');
end

scenarioCRBResults=simCRBResults;
if isRRIS, rrisCRBResults=simCRBResults; else, rrisCRBResults=struct([]); end
save(finalFile,'simCRBResults','rrisCRBResults','scenarioCRBResults', ...
    'params0','crbAlg','criticalAlg', ...
    'searchOptsMultiplicative','searchOptsAdditive', ...
    'gammaSeedSIM','gammaGridByL','formalStepDB','refine','dualState','commBoundary', ...
    'abnormalCRBJumpDB','abnormalEigDropDecades', ...
    'backwardPolishRelTol','enableBackwardPolish', ...
    'enableMonotoneEnvelope','monotoneEnvelopeRelTol', ...
        'enableBackwardReoptimization','backwardRepairTolDB', ...
        'maxBackwardRepairPass','scriptVersion','simCommon','scenarioCRB','-v7.3');
fprintf('\n%s final data saved to:\n%s\n',scenarioCRB,finalFile);
save_transition_diagnostics(diagnosticsFile,transitionDiagnostics, ...
    scriptVersion,dualState,refine,diagnosticOptions);
fprintf('%s transition diagnostics saved to:\n%s\n',scenarioCRB,diagnosticsFile);
fprintf('%s compact per-point diagnostics saved to:\n%s\n', ...
    scenarioCRB,diagnosticOptions.traceTextFile);
fprintf(['After the simulation, run: SIM_print_transition_diagnostics\n', ...
    'to print the complete saved report.\n']);

%% Exact physical CRB versus target and achieved common SINR.
figure('Name',[scenarioCRB ' exact physical CRB versus SINR in dB']);
hold on; grid on; box on;
for idx = 1:numel(simCRBResults)
    tr = simCRBResults(idx).trade;
    ok = arrayfun(@is_valid_trade_record,tr);
    if ~any(ok)
        continue;
    end
    plot([tr(ok).CRBPhysicalDB],[tr(ok).targetGammaDB], ...
        '-o','LineWidth',1.5, ...
        'DisplayName',sprintf('%s L=%d target',scenarioCRB,simCRBResults(idx).L));
    plot([tr(ok).CRBPhysicalDB],[tr(ok).minSINRdB], ...
        '--s','LineWidth',1.0, ...
        'DisplayName',sprintf('%s L=%d achieved',scenarioCRB,simCRBResults(idx).L));
end
xlabel('Sensing CRB (dB)');
ylabel('Common SINR (dB)');
title(sprintf('%s: exact physical CRB--SINR tradeoff',scenarioCRB));
legend('Location','best');

%% Regularized CRB diagnostic versus the same target/achieved SINR points.
% These values are computed only for this figure and are never inserted into
% simCRBResults or saved to the MAT result files.
plotRegAlg = crbAlg;
plotRegAlg.crbInternalRegRel = crbAlg.crbPlotRegRel;
figure('Name',[scenarioCRB ' regularized CRB diagnostic versus SINR in dB']);
hold on; grid on; box on;
for idx = 1:numel(simCRBResults)
    tr = simCRBResults(idx).trade;
    okIdx = find(arrayfun(@is_valid_trade_record,tr));
    if isempty(okIdx)
        continue;
    end
    regCRBPhysicalDB = NaN(1,numel(okIdx));
    targetDB = NaN(1,numel(okIdx));
    achievedDB = NaN(1,numel(okIdx));
    for jj = 1:numel(okIdx)
        ip = okIdx(jj);
        [Jreg,regInfo] = SIM_crb_value_regularized( ...
            simCRBResults(idx).channel,tr(ip).P,tr(ip).Rx,plotRegAlg);
        if regInfo.valid && isfinite(Jreg) && Jreg > 0
            crbScale = tr(ip).sigma_s2/params0.T;
            regCRBPhysicalDB(jj) = ...
                10*log10(max(crbScale*Jreg,realmin));
        end
        targetDB(jj) = tr(ip).targetGammaDB;
        achievedDB(jj) = tr(ip).minSINRdB;
    end
    validReg = isfinite(regCRBPhysicalDB);
    if ~any(validReg)
        continue;
    end
    plot(regCRBPhysicalDB(validReg),targetDB(validReg), ...
        '-o','LineWidth',1.5, ...
        'DisplayName',sprintf('%s L=%d target',scenarioCRB,simCRBResults(idx).L));
    plot(regCRBPhysicalDB(validReg),achievedDB(validReg), ...
        '--s','LineWidth',1.0, ...
        'DisplayName',sprintf('%s L=%d achieved',scenarioCRB,simCRBResults(idx).L));
end
xlabel('Regularized sensing CRB (dB)');
ylabel('Common SINR (dB)');
title(sprintf('%s: regularized CRB diagnostic, regRel = %.0e', ...
    scenarioCRB,crbAlg.crbPlotRegRel));
legend('Location','best');

%% Deferred rank-one recovery (never used by optimization or plotting).
for idx = 1:numel(simCRBResults)
    for ip = 1:numel(simCRBResults(idx).trade)
        recPoint = simCRBResults(idx).trade(ip);
        if ~is_valid_trade_record(recPoint), continue; end
        rec = SIM_recover_CRB_beamformers(params0, ...
            simCRBResults(idx).channel,recPoint.P, ...
            recPoint.QcommSDR,recPoint.RsSDR,crbAlg);
        simCRBResults(idx).trade(ip).rankOneRecovery = rec;
        if rec.valid
            simCRBResults(idx).trade(ip).WRankOne = rec.Wc;
            simCRBResults(idx).trade(ip).R0RankOne = rec.RsRecovered;
        end
    end
end
scenarioCRBResults=simCRBResults;
if isRRIS, rrisCRBResults=simCRBResults; end
save(finalFile,'simCRBResults','rrisCRBResults','scenarioCRBResults', ...
    'params0','crbAlg','criticalAlg', ...
    'searchOptsMultiplicative','searchOptsAdditive', ...
    'gammaSeedSIM','gammaGridByL','formalStepDB','refine','dualState','commBoundary', ...
    'abnormalCRBJumpDB','abnormalEigDropDecades', ...
    'backwardPolishRelTol','enableBackwardPolish', ...
    'enableMonotoneEnvelope','monotoneEnvelopeRelTol', ...
    'enableBackwardReoptimization','backwardRepairTolDB', ...
    'maxBackwardRepairPass','scriptVersion','simCommon','scenarioCRB','-v7.3');
fprintf('Deferred rank-one recovery appended after plotting and saved.\n');
% Automatically create a bounded text report. Detailed traces remain
% available by explicitly running SIM_print_transition_diagnostics('all').
try
    summaryReportFile=fullfile(resultsDir, ...
        [fileTag '_CRB_SINR_diagnostics_report.txt']);
    SIM_print_transition_diagnostics('summary',diagnosticsFile, ...
        finalFile,summaryReportFile);
catch ME
    warning('SIM_CRB_SINR:RecoverySummaryReportFailed', ...
        'Could not write the compact diagnostics report: %s',ME.message);
end

%% ----------------------------- Local functions -------------------------
function [bestCandidate,candidates,severity,reason,transition,carrierInfo, ...
    nextCarrierInit] = ...
    solve_compact_formal_target(params,ch,alg,warmInit,previousRecord, ...
    gammaValue,gammaTarget,refine,seedIndex,useCarrier, ...
    previousCarrierInit) %#ok<INUSD>
% Every target compares independent best-CRB and persistent carrier paths.
% Only complete, exact-CRB/SINR-feasible AO results may enter trade.

candidates = [];
nextCarrierInit=[];
carrierInfo = struct('attempted',false,'selected',false, ...
    'handoffSource','report','status','not-attempted','numAttempts',0, ...
    'blendWeight',NaN,'nextMargin',NaN,'carrierCRB',NaN, ...
    'relMinEigA',NaN,'relMinEigB',NaN, ...
    'targetCurrentMargin',NaN,'targetNextMargin',NaN, ...
    'crbDeltaDB',NaN,'quality','none','usedPreviousCarrier',false, ...
    'buildDetails',struct(),'buildElapsedSeconds',0, ...
    'primaryPath',empty_ao_path_diagnostics('continuation'), ...
    'carrierPath',empty_ao_path_diagnostics('commCarrier'), ...
    'carrierPathSkipReason','not-attempted','additionalCandidates',[]);
carrierInit=[];
if useCarrier
    carrierInfo.attempted=true;
    carrierBuildClock=tic;
    [carrierInit,buildInfo]=SIM_build_comm_carrier(params,ch,alg, ...
        warmInit,previousRecord.targetGammaScalar,gammaValue, ...
        previousCarrierInit);
    carrierInfo.buildElapsedSeconds=toc(carrierBuildClock);
    carrierInfo.buildDetails=buildInfo;
    carrierInfo.status=buildInfo.status;
    carrierInfo.numAttempts=buildInfo.phaseIAttempts;
    carrierInfo.blendWeight=buildInfo.selectedBlendWeight;
    carrierInfo.nextMargin=buildInfo.bestNextMargin;
    carrierInfo.carrierCRB=buildInfo.selectedCarrierCRB;
    carrierInfo.targetCurrentMargin=buildInfo.targetCurrentMargin;
    carrierInfo.targetNextMargin=buildInfo.targetNextMargin;
    carrierInfo.crbDeltaDB=buildInfo.selectedCRBDeltaDB;
    carrierInfo.quality=buildInfo.selectedQuality;
    carrierInfo.usedPreviousCarrier=buildInfo.usedPreviousCarrier;
    if ~isempty(carrierInit)&&buildInfo.successReady
        nextCarrierInit=carrierInit;
        [Pcarrier,~]=build_P(carrierInit.theta0,ch.Omega);
        Rxcarrier=carrierInit.Rs0;
        for kk=1:params.K
            Rxcarrier=Rxcarrier+carrierInit.Qcomm0(:,:,kk);
        end
        [~,carrierCRBInfo]=SIM_crb_value(ch,Pcarrier,Rxcarrier,alg);
        carrierInfo.relMinEigA=carrierCRBInfo.relMinEigA;
        carrierInfo.relMinEigB=carrierCRBInfo.relMinEigB;
    end
end

primaryClock=tic;
[solPrimary,histPrimary] = SIM_run_CRB_SINR_AO_solver( ...
    params,ch,alg,warmInit,gammaTarget);
carrierInfo.primaryPath=summarize_ao_path('continuation', ...
    solPrimary,histPrimary,gammaValue,params.T,toc(primaryClock), ...
    alg.maxOuter);
candidates = append_solution_candidate(candidates, ...
    'continuation',solPrimary,histPrimary,gammaValue);

if useCarrier && ~isempty(nextCarrierInit) && ...
        ~same_complete_init(nextCarrierInit,warmInit)
    carrierClock=tic;
    [solCarrier,histCarrier] = SIM_run_CRB_SINR_AO_solver( ...
        params,ch,alg,nextCarrierInit,gammaTarget);
    carrierInfo.carrierPath=summarize_ao_path('commCarrier', ...
        solCarrier,histCarrier,gammaValue,params.T,toc(carrierClock), ...
        alg.maxOuter);
    carrierInfo.carrierPathSkipReason='';
    candidates = append_solution_candidate(candidates, ...
        'commCarrier',solCarrier,histCarrier,gammaValue);
elseif ~useCarrier
    carrierInfo.carrierPathSkipReason='dual-state-disabled';
elseif isempty(nextCarrierInit)
    carrierInfo.carrierPathSkipReason=['carrier-not-ready:',carrierInfo.status];
else
    carrierInfo.carrierPathSkipReason='same-complete-state-as-bestCRB';
end
if useCarrier
    fprintf(['    commCarrier: %s, margin=%+.3e/target %.3e, ', ...
        'CRBdelta=%+.2f dB (%s), blend=%.3e, inherited=%d\n'], ...
        carrierInfo.status,carrierInfo.nextMargin, ...
        carrierInfo.targetNextMargin,carrierInfo.crbDeltaDB, ...
        carrierInfo.quality,carrierInfo.blendWeight, ...
        carrierInfo.usedPreviousCarrier);
    fprintf(['    paths: bestCRB=%s/%.3f dB/CVX=%s; ', ...
        'carrier=%s/%.3f dB/CVX=%s; time=%.1f/%.1f/%.1fs\n'], ...
        compact_path_status(carrierInfo.primaryPath), ...
        carrierInfo.primaryPath.physicalCRBDB, ...
        carrierInfo.primaryPath.firstCVXStatus, ...
        compact_path_status(carrierInfo.carrierPath), ...
        carrierInfo.carrierPath.physicalCRBDB, ...
        carrierInfo.carrierPath.firstCVXStatus, ...
        carrierInfo.buildElapsedSeconds, ...
        carrierInfo.primaryPath.elapsedSeconds, ...
        carrierInfo.carrierPath.elapsedSeconds);
end

[bestCandidate,found] = select_best_complete_candidate(candidates);
if ~found
    bestCandidate = [];
    severity = 'failure';
    reason = status_text(solPrimary,'status','continuation failed');
    transition = empty_transition();
    return;
end
carrierInfo.selected = strcmp(bestCandidate.label,'commCarrier');
if carrierInfo.selected, carrierInfo.handoffSource='commCarrier'; end
[severity,reason,transition] = classify_solution_transition( ...
    bestCandidate.sol,previousRecord,params,refine);
end

function info=empty_ao_path_diagnostics(label)
info=struct('label',char(label),'attempted',false,'success',false, ...
    'status','not-run','failureStage','not-run', ...
    'physicalCRBDB',NaN,'minSINRDB',NaN,'relativeSINRMargin',NaN, ...
    'firstSDRFeasible',false,'firstSDRStatus','not-run', ...
    'firstCVXStatus','not-run','lastCVXStatus','not-run', ...
    'finalSDRStatus','not-run','cvxFailureCount',0, ...
    'cvxInaccurateCount',0,'sdrRetryCount',0,'outerIterations',0, ...
    'outerIterationBudget',NaN,'hitOuterBudget',false, ...
    'outerStopTriggered',false,'lastRelativeDecrease',NaN, ...
    'bestIteration',NaN,'bestStage','none','initialFeasible',false, ...
    'initialMargin',NaN,'preOuterMarginStart',NaN, ...
    'preOuterMarginEnd',NaN,'wAccepted',0,'phaseAccepted',0, ...
    'innerPhaseSteps',0,'outerMarginSteps',0, ...
    'phaseFeasibilityRejects',0,'phaseObjectiveRejects',0, ...
    'phaseConditionRejects',0,'phaseModelRejects',0, ...
    'phaseLineSearchTrials',0,'maxPhaseGradNorm',NaN, ...
    'lastPhaseGradNorm',NaN,'minRelEigA',NaN,'minRelEigB',NaN, ...
    'elapsedSeconds',0);
end

function info=summarize_ao_path( ...
    label,sol,hist,gammaValue,Tused,elapsed,outerBudget)
if nargin<7,outerBudget=NaN;end
info=empty_ao_path_diagnostics(label);
info.attempted=true;
info.success=is_solution_success(sol);
info.status=status_text(sol,'status','unknown');
info.failureStage=status_text(hist,'failureStage','unknown');
info.firstSDRFeasible=logical_scalar_field(hist,'firstSDRFeasible',false);
info.firstSDRStatus=status_text(hist,'firstSDRStatus','not-run');
info.finalSDRStatus=status_text(hist,'finalSDRStatus','not-run');
info.outerIterations=scalar_field(hist,'outerIterationsCompleted',0);
info.outerIterationBudget=outerBudget;
info.hitOuterBudget=isfinite(outerBudget)&& ...
    info.outerIterations>=outerBudget;
info.outerStopTriggered=logical_scalar_field(hist,'outerStopTriggered',false);
if isfield(hist,'outerRelativeDecrease')
    decreases=hist.outerRelativeDecrease( ...
        isfinite(hist.outerRelativeDecrease));
    if ~isempty(decreases),info.lastRelativeDecrease=decreases(end);end
end
info.bestIteration=scalar_field(sol,'bestIteration',NaN);
info.bestStage=status_text(sol,'bestStage','none');
info.initialFeasible=logical_scalar_field( ...
    hist,'initialIncumbentFeasibleAtTarget',false);
initialSINR=scalar_field(hist,'initialIncumbentMinSINR',NaN);
info.initialMargin=initialSINR/max(gammaValue,realmin)-1;
info.preOuterMarginStart=scalar_field(hist,'preOuterMarginStart',NaN);
info.preOuterMarginEnd=scalar_field(hist,'preOuterMarginEnd',NaN);
info.sdrRetryCount=finite_history_sum(hist,'sdrRetryCount')+ ...
    scalar_field(hist,'finalSDRRetryCount',0);
info.wAccepted=finite_history_sum(hist,'WBlockAccepted');
info.phaseAccepted=finite_history_sum(hist,'phaseBlockAccepted');
info.innerPhaseSteps=finite_history_sum(hist,'phaseInner');
info.outerMarginSteps=finite_history_sum(hist,'outerMarginAcceptedSteps');
info.phaseFeasibilityRejects=finite_history_sum( ...
    hist,'phaseFeasibilityRejects');
info.phaseObjectiveRejects=finite_history_sum(hist,'phaseObjectiveRejects');
info.phaseConditionRejects=finite_history_sum(hist,'phaseConditionRejects');
info.phaseModelRejects=finite_history_sum(hist,'phaseModelRejects');
info.phaseLineSearchTrials=finite_history_sum(hist,'phaseLineSearchTrials');
if isfield(hist,'phaseGradNorm')
    finiteGrad=hist.phaseGradNorm(isfinite(hist.phaseGradNorm));
    if ~isempty(finiteGrad)
        info.maxPhaseGradNorm=max(finiteGrad);
        info.lastPhaseGradNorm=finiteGrad(end);
    end
end
if isfield(hist,'sdrCVXStatus')&&iscell(hist.sdrCVXStatus)
    statuses=hist.sdrCVXStatus(~cellfun(@isempty,hist.sdrCVXStatus));
    if ~isempty(statuses)
        info.firstCVXStatus=char(string(statuses{1}));
        info.lastCVXStatus=char(string(statuses{end}));
        normalized=lower(string(statuses));
        info.cvxFailureCount=sum(contains(normalized,'infeasible')| ...
            contains(normalized,'failed')|contains(normalized,'unbounded'));
        info.cvxInaccurateCount=sum(contains(normalized,'inaccurate'));
    end
end
if info.success
    info.physicalCRBDB=10*log10(max( ...
        sol.sigma_s2/max(Tused,1)*sol.metrics.CRB,realmin));
    minSINR=min(sol.metrics.sinr);
    info.minSINRDB=10*log10(max(minSINR,realmin));
    info.relativeSINRMargin=minSINR/max(gammaValue,realmin)-1;
    info.minRelEigA=sol.metrics.crbInfo.relMinEigA;
    info.minRelEigB=sol.metrics.crbInfo.relMinEigB;
end
info.elapsedSeconds=elapsed;
end

function total=finite_history_sum(hist,fieldName)
total=0;
if ~isstruct(hist)||~isfield(hist,fieldName),return;end
values=double(hist.(fieldName));
values=values(isfinite(values));
if ~isempty(values),total=sum(values);end
end

function text=compact_path_status(info)
if ~isstruct(info)||~isfield(info,'attempted')||~info.attempted
    text='skip';
elseif info.success
    text='ok';
else
    text='fail';
end
end

function summaries=summarize_additional_candidate_paths(candidates,gamma,Tused)
summaries=[];
for ii=1:numel(candidates)
    if ismember(candidates(ii).label,{'continuation','commCarrier'})
        continue;
    end
    item=summarize_ao_path(candidates(ii).label,candidates(ii).sol, ...
        candidates(ii).hist,gamma,Tused,NaN);
    if isempty(summaries),summaries=item;else,summaries(end+1,1)=item;end %#ok<AGROW>
end
end

function [severity,reason,transition] = classify_solution_transition( ...
    sol,previousRecord,params,refine)
transition = empty_transition();
severity = 'none';
reason = '';

if ~is_solution_success(sol)
    severity = 'failure';
    reason = status_text(sol,'status','formal target failed');
    return;
end
if ~is_valid_trade_record(previousRecord)
    return;
end

physicalCRB = (sol.sigma_s2/params.T)*sol.metrics.CRB;
currentDB = 10*log10(max(physicalCRB,realmin));
transition.deltaCRBDB = currentDB-previousRecord.CRBPhysicalDB;
info = sol.metrics.crbInfo;
transition.pointDropA = previousRecord.relMinEigA / ...
    max(info.relMinEigA,realmin);
transition.pointDropB = previousRecord.relMinEigB / ...
    max(info.relMinEigB,realmin);
transition.conditionWarning = logical_scalar_field( ...
    info,'conditionWarning',false);

hardJump = transition.deltaCRBDB > refine.hardCRBJumpDB || ...
    transition.pointDropA > refine.hardEigDropFactor || ...
    transition.pointDropB > refine.hardEigDropFactor;
softJump = transition.deltaCRBDB > refine.softCRBJumpDB || ...
    transition.pointDropA > refine.softEigDropFactor || ...
    transition.pointDropB > refine.softEigDropFactor || ...
    transition.conditionWarning;

if hardJump
    severity = 'hard';
elseif softJump
    severity = 'soft';
else
    return;
end

reason = sprintf(['%s transition: dCRB=%+.2f dB, ', ...
    'pointDropA/B=%.2e/%.2e, conditionWarning=%d'], ...
    severity,transition.deltaCRBDB,transition.pointDropA, ...
    transition.pointDropB,transition.conditionWarning);
end

function transition = empty_transition()
transition.deltaCRBDB = NaN;
transition.pointDropA = NaN;
transition.pointDropB = NaN;
transition.conditionWarning = false;
end

function [bestCandidate,candidates] = solve_branch_switch_target( ...
    params,ch,alg,initBase,warmInit,previousLayerTrade, ...
    gammaValue,gammaTarget,recoveryOpts,refine,seedIndex,depth)
% Branch switching is used only after the local-refinement budget is
% exhausted. The compact builder returns at most two candidates: a random-
% phase warm start and one structured cross-layer (or zero-phase) start.

previousRef = nearest_valid_trade(previousLayerTrade,gammaValue);
if isempty(previousRef)
    previousSol = [];
else
    previousSol = trade_record_to_solution(previousRef);
end

opts = recoveryOpts;
opts.profile = 'compact';
opts.maxCandidates = 2;
opts.randomizeW = false;
opts.seed = recoveryOpts.seed + 10000*ch.L + 100*seedIndex + depth;

baseForRecovery = warmInit;
if ~isfield(baseForRecovery,'theta0') || isempty(baseForRecovery.theta0)
    baseForRecovery = initBase;
end
[initList,labels] = SIM_build_cross_layer_inits( ...
    params,ch,baseForRecovery,previousSol,opts);

candidates = [];
for ii = 1:numel(initList)
    [solTry,histTry] = SIM_run_CRB_SINR_AO_solver( ...
        params,ch,alg,initList{ii},gammaTarget);
    candidates = append_solution_candidate(candidates, ...
        ['branch-',labels{ii}],solTry,histTry,gammaValue);
end
[bestCandidate,found] = select_best_complete_candidate(candidates);
if ~found
    bestCandidate = [];
end
end

function combined = concatenate_solution_candidates(first,second)
if isempty(first)
    combined = second;
    return;
end
if isempty(second)
    combined = first;
    return;
end
combined = first;
for ii = 1:numel(second)
    combined(end+1,1) = orderfields(second(ii),combined(1)); %#ok<AGROW>
end
end

function task = make_refinement_task( ...
    gamma,depth,isInsertedPoint,parentLowGamma,parentHighGamma, ...
    isForcedPoint)
task.gamma = gamma;
task.depth = depth;
task.isInsertedPoint = logical(isInsertedPoint);
task.parentLowGamma = parentLowGamma;
task.parentHighGamma = parentHighGamma;
task.isForcedPoint = logical(isForcedPoint);
end

function suffix = inserted_suffix(isInsertedPoint)
if isInsertedPoint
    suffix = ' [inserted]';
else
    suffix = '';
end
end

function suffix = forced_suffix(isForcedPoint)
if isForcedPoint
    suffix = ' [forced-window]';
else
    suffix = '';
end
end

function record = last_valid_trade_record(trade)
record = empty_trade_record();
for ii = numel(trade):-1:1
    if is_valid_trade_record(trade(ii))
        record = trade(ii);
        return;
    end
end
error('SIM_CRB_SINR:NoValidPreviousTrade', ...
    'No valid previous formal point is available for continuation.');
end

function candidates = append_solution_candidate( ...
    candidates,label,sol,hist,trialGamma)
entry.label = char(string(label));
entry.sol = sol;
entry.hist = hist;
entry.trialGamma = trialGamma;
entry.selectable = true;
entry.success = is_solution_success(sol);
entry.CRB = candidate_crb(sol);
entry.minSINR = NaN;
entry.relMinEigA = NaN;
entry.relMinEigB = NaN;
entry.status = status_text(sol,'status');
if entry.success
    entry.minSINR = min(sol.metrics.sinr);
    entry.relMinEigA = sol.metrics.crbInfo.relMinEigA;
    entry.relMinEigB = sol.metrics.crbInfo.relMinEigB;
end
if isempty(candidates)
    candidates = entry;
else
    candidates(end+1,1) = orderfields(entry,candidates(1)); %#ok<AGROW>
end
end

function [best,found] = select_best_complete_candidate(candidates)
best = [];
found = false;
bestCRB = Inf;
for ii = 1:numel(candidates)
    selectable = true;
    if isfield(candidates(ii),'selectable')
        selectable = logical(candidates(ii).selectable);
    end
    if ~selectable || ~candidates(ii).success || ...
            ~isfinite(candidates(ii).CRB) || candidates(ii).CRB <= 0
        continue;
    end
    if candidates(ii).CRB < bestCRB
        bestCRB = candidates(ii).CRB;
        best = candidates(ii);
        found = true;
    end
end
end

function summary = candidates_to_summary(candidates)
summary = [];
for ii = 1:numel(candidates)
    entry.label = candidates(ii).label;
    entry.trialGamma = candidates(ii).trialGamma;
    entry.success = candidates(ii).success;
    entry.selectable = candidates(ii).selectable;
    entry.CRB = candidates(ii).CRB;
    entry.minSINR = candidates(ii).minSINR;
    entry.relMinEigA = candidates(ii).relMinEigA;
    entry.relMinEigB = candidates(ii).relMinEigB;
    entry.status = candidates(ii).status;
    if isempty(summary)
        summary = entry;
    else
        summary(end+1,1) = orderfields(entry,summary(1)); %#ok<AGROW>
    end
end
end

function print_candidate_summary(label,sol,params)
if ~is_solution_success(sol)
    fprintf('    %-42s failed: %s\n',label,status_text(sol,'status'));
    return;
end
physicalCRB = (sol.sigma_s2/params.T)*sol.metrics.CRB;
info = sol.metrics.crbInfo;
fprintf(['    %-42s CRB=%8.3f dB, SINR=%7.3f dB, ', ...
    'relEig=%.1e/%.1e\n'],label, ...
    10*log10(max(physicalCRB,realmin)), ...
    10*log10(max(min(sol.metrics.sinr),realmin)), ...
    info.relMinEigA,info.relMinEigB);
end

function value = candidate_crb(sol)
value = Inf;
if is_solution_success(sol)
    value = sol.metrics.CRB;
end
end

function tf = is_solution_success(sol)
tf = isstruct(sol) && isfield(sol,'success') && isscalar(sol.success) ...
    && logical(sol.success) && isfield(sol,'metrics') ...
    && isfield(sol.metrics,'CRB') && isscalar(sol.metrics.CRB) ...
    && isfinite(sol.metrics.CRB) && sol.metrics.CRB > 0 ...
    && isfield(sol.metrics,'crbValid') && logical(sol.metrics.crbValid) ...
    && isfield(sol.metrics,'sinrFeasible') && logical(sol.metrics.sinrFeasible);
end

function init = solution_to_init(baseInit,sol,sigma_c2,sigma_s2)
init = baseInit;
if isfield(sol,'W') && ~isempty(sol.W)
    init.W0 = sol.W;
end
if isfield(sol,'QcommSDR') && ~isempty(sol.QcommSDR), init.Qcomm0=sol.QcommSDR; end
if isfield(sol,'RsSDR') && ~isempty(sol.RsSDR), init.Rs0=sol.RsSDR; end
if isfield(sol,'theta') && ~isempty(sol.theta)
    init.theta0 = sol.theta;
end
if isfield(sol,'S') && ~isempty(sol.S)
    init.S = sol.S;
end
init.sigma_c2 = sigma_c2;
init.sigma_s2 = sigma_s2;
end

function init = trade_record_to_init(baseInit,record,sigma_c2,sigma_s2)
init = baseInit;
if isfield(record,'W') && ~isempty(record.W)
    init.W0 = record.W;
end
if isfield(record,'QcommSDR') && ~isempty(record.QcommSDR), init.Qcomm0=record.QcommSDR; end
if isfield(record,'RsSDR') && ~isempty(record.RsSDR), init.Rs0=record.RsSDR; end
if isfield(record,'theta') && ~isempty(record.theta)
    init.theta0 = record.theta;
end
if isfield(record,'S') && ~isempty(record.S)
    init.S = record.S;
end
init.sigma_c2 = sigma_c2;
init.sigma_s2 = sigma_s2;
end

function sol = trade_record_to_solution(record)
sol.success = is_valid_trade_record(record);
sol.status = record.status;
sol.W = record.W;
sol.Wc = record.Wc;
sol.Wr = record.Wr;
sol.theta = record.theta;
sol.P = record.P;
sol.Rx = record.Rx;
sol.QcommSDR = record.QcommSDR;
sol.RsSDR = record.RsSDR;
sol.RsRecovered = record.RsRecovered;
sol.S = record.S;
sol.sigma_c2 = record.sigma_c2;
sol.sigma_s2 = record.sigma_s2;
sol.metrics.CRB = record.CRB;
sol.metrics.sinr = record.sinr;
sol.metrics.MIc = record.MIc;
sol.metrics.MIs = record.MIs;
sol.metrics.rate = record.rate;
sol.metrics.nmmse = record.nmmse;
sol.metrics.crbInfo.rcondA = record.rcondA;
sol.metrics.crbInfo.rcondB = record.rcondB;
sol.metrics.crbInfo.relMinEigA = record.relMinEigA;
sol.metrics.crbInfo.relMinEigB = record.relMinEigB;
end

function [list,labels] = append_unique_init(list,labels,candidate,label)
for ii = 1:numel(list)
    if same_phase_init(list{ii},candidate)
        return;
    end
end
list{end+1,1} = candidate; %#ok<AGROW>
labels{end+1,1} = label; %#ok<AGROW>
end

function tf = same_phase_init(a,b)
tf = false;
if ~isstruct(a) || ~isstruct(b) || ...
        ~isfield(a,'theta0') || ~isfield(b,'theta0') || ...
        isempty(a.theta0) || isempty(b.theta0) || ...
        ~isequal(size(a.theta0),size(b.theta0))
    return;
end
tf = norm(exp(1j*a.theta0(:))-exp(1j*b.theta0(:))) <= ...
    1e-12*sqrt(numel(a.theta0));
end

function tf=same_complete_init(a,b)
tf=same_phase_init(a,b);
if ~tf||~isfield(a,'Qcomm0')||~isfield(b,'Qcomm0')|| ...
        ~isfield(a,'Rs0')||~isfield(b,'Rs0')|| ...
        ~isequal(size(a.Qcomm0),size(b.Qcomm0))|| ...
        ~isequal(size(a.Rs0),size(b.Rs0))
    tf=false;
    return;
end
qScale=max([norm(a.Qcomm0(:)),norm(b.Qcomm0(:)),realmin]);
rScale=max([norm(a.Rs0(:)),norm(b.Rs0(:)),realmin]);
tf=norm(a.Qcomm0(:)-b.Qcomm0(:))<=1e-10*qScale&& ...
    norm(a.Rs0(:)-b.Rs0(:))<=1e-10*rScale;
end

function record = empty_trade_record()
record.success = false;
record.status = '';
record.targetGamma = [];
record.targetGammaScalar = NaN;
record.targetGammaDB = NaN;
record.targetRatePerUser = [];
record.targetRateSum = NaN;
record.gammaSeed = NaN;
record.gammaMaxFound = NaN;
record.CRB = NaN;
record.CRBdB = NaN;
record.CRBObjective = NaN;
record.CRBObjectiveDB = NaN;
record.CRBScale = NaN;
record.CRBPhysical = NaN;
record.CRBPhysicalDB = NaN;
record.MIc = NaN;
record.MIs = NaN;
record.rate = NaN;
record.nmmse = NaN;
record.sinr = [];
record.minSINR = NaN;
record.minSINRdB = NaN;
record.maxScaledConstraintViolation = NaN;
record.rcondA = NaN;
record.rcondB = NaN;
record.minEigA = NaN;
record.maxEigA = NaN;
record.relMinEigA = NaN;
record.minEigB = NaN;
record.maxEigB = NaN;
record.relMinEigB = NaN;
record.traceAinv = NaN;
record.traceBinv = NaN;
record.logCRB = NaN;
record.conditionWarning = false;
record.bestIteration = NaN;
record.bestStage = 'none';
record.failureStage = 'not-reported';
record.firstSDRFeasible = false;
record.W = [];
record.Wc = [];
record.Wr = [];
record.theta = [];
record.P = [];
record.Rx = [];
record.QcommSDR = [];
record.RsSDR = [];
record.WUnrecovered = [];
record.R0Unrecovered = [];
record.RsRecovered = [];
record.rankOneRecovery = [];
record.WRankOne = [];
record.R0RankOne = [];
record.S = [];
record.sigma_c2 = NaN;
record.sigma_s2 = NaN;
record.hist = [];
record.selectedInitialization = 'none';
record.initializationTrials = [];
record.abnormalRecoveryTriggered = false;
record.abnormalRecoveryReason = '';
record.isInsertedPoint = false;
record.isForcedPoint = false;
record.useCriticalBudget = false;
record.maxOuterUsed = NaN;
record.maxPhaseInnerUsed = NaN;
record.refinementDepth = 0;
record.parentLowGamma = NaN;
record.parentHighGamma = NaN;
record.transitionSeverity = 'none';
record.deltaCRBDBFromPrevious = NaN;
record.pointDropA = NaN;
record.pointDropB = NaN;
record.branchWarning = false;
record.branchRejected = false;
record.branchProtectionExhausted = false;
record.branchSwitchTried = false;
record.phaseGuardTriggered = false;
record.phaseGuardRejectCount = 0;
record.phaseAbsoluteGuardRejectCount = 0;
record.phaseRelativeGuardRejectCount = 0;
record.phaseGuardMaxDropA = 1;
record.phaseGuardMaxDropB = 1;
record.phaseGuardMinRelEigA = NaN;
record.phaseGuardMinRelEigB = NaN;
record.recoveryValid = false;
record.recoveryStatus = 'not-reported';
record.recoveryCovarianceGapRel = NaN;
record.recoveryPowerGapRel = NaN;
record.recoveryDesiredGapRel = NaN;
record.recoverySelfLeakRel = NaN;
record.recoveryMinResidualRelEig = NaN;
record.solutionSource = 'direct-solve';
record.sourceTargetGamma = NaN;
record.incumbentSourcePoint = NaN;
record.polishSourcePoint = NaN;
record.backwardReoptimized = false;
record.backwardSourcePoint = NaN;
end

function record = make_trade_record(sol,hist,gammaValue,gammaTarget, ...
    gammaSeed,gammaMaxFound,params,selectedInitLabel,trialSummary, ...
    recoveryTriggered,recoveryReason)
record = empty_trade_record();
record.targetGamma = gammaTarget;
record.targetGammaScalar = gammaValue;
record.targetGammaDB = 10*log10(max(gammaValue,realmin));
record.targetRatePerUser = log2(1+gammaTarget);
record.targetRateSum = sum(record.targetRatePerUser);
record.gammaSeed = gammaSeed;
record.gammaMaxFound = gammaMaxFound;
record.selectedInitialization = selectedInitLabel;
record.initializationTrials = trialSummary;
record.abnormalRecoveryTriggered = recoveryTriggered;
record.abnormalRecoveryReason = recoveryReason;
record.status = status_text(sol,'status');
record.success = is_solution_success(sol);
record.phaseGuardTriggered = logical_scalar_field( ...
    sol,'phaseGuardTriggered',false);
record.phaseGuardRejectCount = scalar_field( ...
    sol,'phaseGuardRejectCount',0);
record.phaseAbsoluteGuardRejectCount = scalar_field( ...
    sol,'phaseAbsoluteGuardRejectCount',0);
record.phaseRelativeGuardRejectCount = scalar_field( ...
    sol,'phaseRelativeGuardRejectCount',0);
record.phaseGuardMaxDropA = scalar_field(sol,'phaseGuardMaxDropA',1);
record.phaseGuardMaxDropB = scalar_field(sol,'phaseGuardMaxDropB',1);
record.phaseGuardMinRelEigA = scalar_field( ...
    sol,'phaseGuardMinRelEigA',NaN);
record.phaseGuardMinRelEigB = scalar_field( ...
    sol,'phaseGuardMinRelEigB',NaN);
record.recoveryValid = logical_scalar_field(sol,'recoveryValid',false);
record.recoveryStatus = status_text(sol,'recoveryStatus','not-reported');
record.recoveryCovarianceGapRel = scalar_field( ...
    sol,'recoveryCovarianceGapRel',NaN);
record.recoveryPowerGapRel = scalar_field(sol,'recoveryPowerGapRel',NaN);
record.recoveryDesiredGapRel = scalar_field( ...
    sol,'recoveryDesiredGapRel',NaN);
record.recoverySelfLeakRel = scalar_field(sol,'recoverySelfLeakRel',NaN);
record.recoveryMinResidualRelEig = scalar_field( ...
    sol,'recoveryMinResidualRelEig',NaN);

if ~record.success
    record.failureStage = status_text(sol,'failureStage');
    if isfield(sol,'firstSDRFeasible') && isscalar(sol.firstSDRFeasible)
        record.firstSDRFeasible = logical(sol.firstSDRFeasible);
    end
    return;
end

crbObjective = sol.metrics.CRB;
crbScale = sol.sigma_s2/params.T;
crbPhysical = crbScale*crbObjective;
minSINR = min(sol.metrics.sinr);
info = sol.metrics.crbInfo;

record.CRB = crbObjective;
record.CRBdB = 10*log10(max(crbObjective,realmin));
record.CRBObjective = crbObjective;
record.CRBObjectiveDB = record.CRBdB;
record.CRBScale = crbScale;
record.CRBPhysical = crbPhysical;
record.CRBPhysicalDB = 10*log10(max(crbPhysical,realmin));
record.MIc = sol.metrics.MIc;
record.MIs = sol.metrics.MIs;
record.rate = sol.metrics.rate;
record.nmmse = sol.metrics.nmmse;
record.sinr = sol.metrics.sinr;
record.minSINR = minSINR;
record.minSINRdB = 10*log10(max(minSINR,realmin));
record.maxScaledConstraintViolation = ...
    sol.metrics.maxScaledConstraintViolation;
record.rcondA = info.rcondA;
record.rcondB = info.rcondB;
record.minEigA = info.minEigA;
record.maxEigA = info.maxEigA;
record.relMinEigA = info.relMinEigA;
record.minEigB = info.minEigB;
record.maxEigB = info.maxEigB;
record.relMinEigB = info.relMinEigB;
record.traceAinv = info.traceAinv;
record.traceBinv = info.traceBinv;
record.logCRB = info.logCRB;
record.conditionWarning = logical(info.conditionWarning);
record.bestIteration = scalar_field(sol,'bestIteration',NaN);
record.bestStage = status_text(sol,'bestStage','none');
record.failureStage = status_text(sol,'failureStage','none');
record.firstSDRFeasible = logical_scalar_field( ...
    sol,'firstSDRFeasible',false);
record.W = sol.W;
record.Wc = sol.Wc;
record.Wr = sol.Wr;
record.theta = sol.theta;
record.P = sol.P;
record.Rx = sol.Rx;
record.QcommSDR = sol.QcommSDR;
record.RsSDR = sol.RsSDR;
record.WUnrecovered = sol.QcommSDR;
record.R0Unrecovered = sol.RsSDR;
record.RsRecovered = sol.RsRecovered;
record.S = sol.S;
record.sigma_c2 = sol.sigma_c2;
record.sigma_s2 = sol.sigma_s2;
record.hist = hist;
end

function sol = point_failure_solution(gammaTarget)
sol.success = false;
sol.status = 'All formal-point initializations failed.';
sol.failureStage = 'all-initializations-failed';
sol.firstSDRFeasible = false;
sol.gammaTarget = gammaTarget;
end

function trade = refresh_trade_transition_metadata(trade,refine)
previous = [];
for ii = 1:numel(trade)
    if ~is_valid_trade_record(trade(ii))
        continue;
    end
    if isempty(previous)
        trade(ii).transitionSeverity = 'none';
        trade(ii).deltaCRBDBFromPrevious = NaN;
        trade(ii).pointDropA = NaN;
        trade(ii).pointDropB = NaN;
        trade(ii).branchWarning = false;
    else
        deltaDB = trade(ii).CRBPhysicalDB-previous.CRBPhysicalDB;
        dropA = previous.relMinEigA/max(trade(ii).relMinEigA,realmin);
        dropB = previous.relMinEigB/max(trade(ii).relMinEigB,realmin);
        hard = deltaDB > refine.hardCRBJumpDB || ...
            dropA > refine.hardEigDropFactor || ...
            dropB > refine.hardEigDropFactor;
        soft = deltaDB > refine.softCRBJumpDB || ...
            dropA > refine.softEigDropFactor || ...
            dropB > refine.softEigDropFactor || ...
            trade(ii).conditionWarning;
        if hard
            severity = 'hard';
        elseif soft
            severity = 'soft';
        else
            severity = 'none';
        end
        trade(ii).transitionSeverity = severity;
        trade(ii).deltaCRBDBFromPrevious = deltaDB;
        trade(ii).pointDropA = dropA;
        trade(ii).pointDropB = dropB;
        trade(ii).branchWarning = ~strcmp(severity,'none');
    end
    previous = trade(ii);
end
end

function out = copy_refinement_metadata(out,source)
if ~isstruct(out) || ~isstruct(source)
    return;
end
fields = {'isInsertedPoint','isForcedPoint','refinementDepth', ...
    'parentLowGamma','parentHighGamma', ...
    'branchProtectionExhausted','branchSwitchTried'};
for ii = 1:numel(fields)
    name = fields{ii};
    if isfield(source,name)
        out.(name) = source.(name);
    end
end
end

function [count,maxGainDB] = diagnose_backward_monotonicity(trade,relTol)
%DIAGNOSE_BACKWARD_MONOTONICITY Report without replacing forward records.
count = 0;
maxGainDB = 0;
for ii = 1:numel(trade)-1
    low = trade(ii);
    high = trade(ii+1);
    if ~is_valid_trade_record(low) || ~is_valid_trade_record(high)
        continue;
    end
    if low.CRBPhysical > high.CRBPhysical*(1+relTol)
        count = count + 1;
        gainDB = low.CRBPhysicalDB-high.CRBPhysicalDB;
        maxGainDB = max(maxGainDB,gainDB);
    end
end
end

function trade = backward_polish_complete_tradeoff( ...
    trade,params,ch,alg,criticalAlg,L,refine, ...
    baseInit,sigma_c2,sigma_s2,gammaSeed,gammaMaxFound,relTol)

for ip = numel(trade)-1:-1:1
    high = trade(ip+1);
    low = trade(ip);
    if ~is_valid_trade_record(high)
        continue;
    end
    if is_valid_trade_record(low) && ...
            low.CRBPhysical <= high.CRBPhysical*(1+relTol)
        continue;
    end

    fprintf('\n  [backward repair %02d <- %02d] target %.3f dB\n', ...
        ip,ip+1,low.targetGammaDB);

    gammaTarget = low.targetGammaScalar*ones(params.K,1);
    candidates = [];

    if is_in_forced_window( ...
            L,low.targetGammaScalar,high.targetGammaScalar,refine)
        algPolish = criticalAlg;
    else
        algPolish = alg;
    end

    if is_valid_trade_record(low)
        candidates = append_record_candidate(candidates, ...
            'existing-lower-point',low);
    end

    polishInit = trade_record_to_init(baseInit,high,sigma_c2,sigma_s2);
    [solPolish,histPolish] = SIM_run_CRB_SINR_AO_solver( ...
        params,ch,algPolish,polishInit,gammaTarget);
    if is_solution_success(solPolish)
        summary = [];
        polished = make_trade_record(solPolish,histPolish, ...
            low.targetGammaScalar,gammaTarget,gammaSeed,gammaMaxFound, ...
            params,sprintf('backward-polish-from-point-%d',ip+1), ...
            summary,false,'');
        polished = copy_refinement_metadata(polished,low);
        polished.useCriticalBudget = ...
            algPolish.maxOuter == criticalAlg.maxOuter && ...
            algPolish.maxPhaseInner == criticalAlg.maxPhaseInner;
        polished.maxOuterUsed = algPolish.maxOuter;
        polished.maxPhaseInnerUsed = algPolish.maxPhaseInner;
        polished.solutionSource = 'backward-reoptimization';
        polished.polishSourcePoint = ip+1;
        candidates = append_record_candidate(candidates, ...
            'backward-reoptimization',polished);
    end

    incumbent = retarget_complete_record(high,low.targetGammaScalar, ...
        gammaTarget,params,ch,algPolish,ip+1);
    incumbent = copy_refinement_metadata(incumbent,low);
    if is_valid_trade_record(incumbent)
        incumbent.useCriticalBudget = high.useCriticalBudget;
        incumbent.maxOuterUsed = high.maxOuterUsed;
        incumbent.maxPhaseInnerUsed = high.maxPhaseInnerUsed;
    end
    if is_valid_trade_record(incumbent)
        candidates = append_record_candidate(candidates, ...
            'higher-target-complete-incumbent',incumbent);
    end

    [bestRecord,found] = select_best_record_candidate(candidates);
    if found && (~is_valid_trade_record(low) || ...
            bestRecord.CRBPhysical < low.CRBPhysical*(1-relTol))
        trade(ip) = bestRecord;
        fprintf(['    retained %s, CRB %.3f dB, achieved SINR %.3f dB\n'], ...
            bestRecord.solutionSource,bestRecord.CRBPhysicalDB, ...
            bestRecord.minSINRdB);
    else
        fprintf('    no complete-solution improvement retained.\n');
    end
end
end

function record = retarget_complete_record( ...
    high,gammaValue,gammaTarget,params,ch,alg,sourcePoint)
record = empty_trade_record();
if ~is_valid_trade_record(high)
    return;
end
% The channel structure built by SIM_build_channels does not store the
% externally calibrated noise powers. During the AO solver they are copied
% from initState into a local channel structure, but this retargeting path
% evaluates a saved complete solution directly and therefore must restore
% the noise values saved with that same solution.
chEval = ch;
if isfield(high,'sigma_c2') && isscalar(high.sigma_c2) && ...
        isfinite(high.sigma_c2) && high.sigma_c2 > 0
    chEval.sigma_c2 = high.sigma_c2;
else
    error('SIM_CRB_SINR:MissingIncumbentSigmaC2', ...
        ['The saved complete incumbent does not contain a valid ', ...
         'communication-noise variance sigma_c2.']);
end
if isfield(high,'sigma_s2') && isscalar(high.sigma_s2) && ...
        isfinite(high.sigma_s2) && high.sigma_s2 > 0
    chEval.sigma_s2 = high.sigma_s2;
else
    error('SIM_CRB_SINR:MissingIncumbentSigmaS2', ...
        ['The saved complete incumbent does not contain a valid ', ...
         'sensing-noise variance sigma_s2.']);
end

metrics = SIM_evaluate_CRB_SINR_metrics( ...
    params,chEval,high.P,high.QcommSDR,high.RsSDR,high.S,high.Rx,gammaTarget,alg);
if ~metrics.crbValid || ~metrics.sinrFeasible
    return;
end
record = high;
record.targetGamma = gammaTarget;
record.targetGammaScalar = gammaValue;
record.targetGammaDB = 10*log10(max(gammaValue,realmin));
record.targetRatePerUser = log2(1+gammaTarget);
record.targetRateSum = sum(record.targetRatePerUser);
record.CRB = metrics.CRB;
record.CRBdB = 10*log10(max(metrics.CRB,realmin));
record.CRBObjective = metrics.CRB;
record.CRBObjectiveDB = record.CRBdB;
record.CRBPhysical = record.CRBScale*metrics.CRB;
record.CRBPhysicalDB = 10*log10(max(record.CRBPhysical,realmin));
record.MIc = metrics.MIc;
record.MIs = metrics.MIs;
record.rate = metrics.rate;
record.nmmse = metrics.nmmse;
record.sinr = metrics.sinr;
record.minSINR = min(metrics.sinr);
record.minSINRdB = 10*log10(max(record.minSINR,realmin));
record.maxScaledConstraintViolation = metrics.maxScaledConstraintViolation;
record.rcondA = metrics.crbInfo.rcondA;
record.rcondB = metrics.crbInfo.rcondB;
record.minEigA = metrics.crbInfo.minEigA;
record.maxEigA = metrics.crbInfo.maxEigA;
record.relMinEigA = metrics.crbInfo.relMinEigA;
record.minEigB = metrics.crbInfo.minEigB;
record.maxEigB = metrics.crbInfo.maxEigB;
record.relMinEigB = metrics.crbInfo.relMinEigB;
record.traceAinv = metrics.crbInfo.traceAinv;
record.traceBinv = metrics.crbInfo.traceBinv;
record.logCRB = metrics.crbInfo.logCRB;
record.conditionWarning = metrics.crbInfo.conditionWarning;
record.selectedInitialization = sprintf( ...
    'complete-incumbent-from-point-%d',sourcePoint);
record.initializationTrials = [];
record.solutionSource = 'higher-target-complete-incumbent';
record.sourceTargetGamma = high.targetGammaScalar;
record.incumbentSourcePoint = sourcePoint;
record.polishSourcePoint = sourcePoint;
record.status = sprintf([ ...
    'Complete W/P solution from point %d re-evaluated at the lower target.'], ...
    sourcePoint);
record.success = true;
end

function candidates = append_record_candidate(candidates,label,record)
entry.label = label;
entry.record = record;
entry.CRBPhysical = record.CRBPhysical;
entry.success = is_valid_trade_record(record);
if isempty(candidates)
    candidates = entry;
else
    candidates(end+1,1) = orderfields(entry,candidates(1)); %#ok<AGROW>
end
end

function [record,found] = select_best_record_candidate(candidates)
record = empty_trade_record();
found = false;
best = Inf;
for ii = 1:numel(candidates)
    if candidates(ii).success && isfinite(candidates(ii).CRBPhysical) ...
            && candidates(ii).CRBPhysical < best
        best = candidates(ii).CRBPhysical;
        record = candidates(ii).record;
        found = true;
    end
end
end

function ref = nearest_valid_trade(trade,gammaValue)
ref = [];
if isempty(trade)
    return;
end
bestDistance = Inf;
for ii = 1:numel(trade)
    if ~is_valid_trade_record(trade(ii))
        continue;
    end
    distance = abs(10*log10(trade(ii).targetGammaScalar) ...
        - 10*log10(gammaValue));
    if distance < bestDistance
        bestDistance = distance;
        ref = trade(ii);
    end
end
end

function tf = is_valid_trade_record(record)
tf = isstruct(record) && isfield(record,'success') ...
    && isscalar(record.success) && logical(record.success) ...
    && isfield(record,'CRBPhysical') && isscalar(record.CRBPhysical) ...
    && isfinite(record.CRBPhysical) && record.CRBPhysical > 0 ...
    && isfield(record,'QcommSDR') && ~isempty(record.QcommSDR) ...
    && isfield(record,'RsSDR') && ~isempty(record.RsSDR) ...
    && isfield(record,'theta') && ~isempty(record.theta) ...
    && isfield(record,'P') && ~isempty(record.P) ...
    && isfield(record,'Rx') && ~isempty(record.Rx);
end

function tf = is_in_forced_window(L,gammaLow,gammaHigh,refine)
tf = false;
if ~isfield(refine,'forceWindowByL') || ...
        L < 1 || L > numel(refine.forceWindowByL)
    return;
end
windowDB = refine.forceWindowByL{L};
if isempty(windowDB)
    return;
end
if ~isnumeric(windowDB) || numel(windowDB) ~= 2 || ...
        any(~isfinite(windowDB)) || windowDB(2) <= windowDB(1)
    error('SIM_CRB_SINR:InvalidForcedWindow', ...
        'refine.forceWindowByL{%d} must be [lowDB, highDB].',L);
end
lowDB = 10*log10(max(gammaLow,realmin));
highDB = 10*log10(max(gammaHigh,realmin));
if highDB < lowDB
    tmp = lowDB;
    lowDB = highDB;
    highDB = tmp;
end
% Half-open overlap: [lowDB,highDB] intersects [windowLow,windowHigh].
% This includes 2--3,...,8--9 dB for the configured [2,9] dB window.
tf = highDB > windowDB(1) && lowDB < windowDB(2);
end

function tf = is_in_dual_state_window(L,gammaLow,gammaHigh,dualState)
tf = false;
if ~isstruct(dualState) || ~isfield(dualState,'enable') || ...
        ~logical(dualState.enable) || gammaHigh <= gammaLow
    return;
end
if ~isfield(dualState,'windowByL') || L > numel(dualState.windowByL)
    return;
end
windowDB = dualState.windowByL{L};
if isempty(windowDB) || numel(windowDB) ~= 2
    return;
end
lowDB = 10*log10(max(gammaLow,realmin));
highDB = 10*log10(max(gammaHigh,realmin));
tf = highDB >= min(windowDB) && lowDB <= max(windowDB);
end

function gammaGrid = build_uniform_db_grid(gammaStart,gammaMax,stepDB)
validateattributes(gammaStart,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(gammaMax,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(stepDB,{'numeric'}, ...
    {'real','finite','scalar','positive'});
if gammaMax <= gammaStart*(1+1e-12)
    gammaGrid = gammaMax;
    return;
end
startDB = 10*log10(gammaStart);
maxDB = 10*log10(gammaMax);
gridDB = startDB:stepDB:maxDB;
if isempty(gridDB)
    gridDB = [startDB,maxDB];
elseif gridDB(end) < maxDB-1e-10
    gridDB(end+1) = maxDB;
else
    gridDB(end) = maxDB;
end
gammaGrid = 10.^(gridDB(:)/10);
gammaGrid(1) = gammaStart;
gammaGrid(end) = gammaMax;
end

function entry = make_SIM_result_entry( ...
    L,mode,ch,initBase,gammaSeed,gammaRequestedRangeStart, ...
    gammaMaxFound,gammaGrid,formalStepDB,searchInfo,rangeTrials, ...
    selectedRangeInitLabel,trade,noiseInfo,sigma_c2,sigma_s2,T)
if isfield(ch,'type'), entry.scenario=ch.type; else, entry.scenario='SIM'; end
entry.L = L;
entry.mode = mode;
entry.channel = ch;
entry.initialW0 = initBase.W0;
entry.initialTheta0 = initBase.theta0;
entry.gammaSeed = gammaSeed;
entry.gammaRequestedRangeStart = gammaRequestedRangeStart;
entry.gammaMaxFound = gammaMaxFound;
entry.gammaGrid = gammaGrid;
entry.formalStepDB = formalStepDB;
entry.tradeRelativeRatio = 10^(formalStepDB/10);
entry.tradeAbsoluteStep = NaN;
entry.gammaSearchInfo = searchInfo;
entry.rangeTrials = rangeTrials;
entry.selectedRangeInitLabel = selectedRangeInitLabel;
entry.trade = trade;
entry.commCarrier=cell(numel(trade),1);
entry.commCarrierDebug=cell(numel(trade),1);
% These fields are overwritten with the actual raw/enveloped curves after
% the forward layer completes. Defaults keep checkpoint structures uniform.
entry.tradeForwardRaw = trade;
% Keep the top-level result-entry schema fixed from the very first
% checkpoint assignment. These fields are populated with the actual
% backward-reoptimized data after the layer finishes.
entry.tradeBackwardReoptimized = trade;
entry.backwardRepairInfo = struct('enabled',false);
entry.tradeMonotone = trade;
entry.monotoneInfo = struct('enabled',false,'numReplacements',0, ...
    'maxImprovementDB',0,'rawViolationCount',0, ...
    'finalViolationCount',0);
entry.noiseInfo = noiseInfo;
entry.sigma_c2_used = sigma_c2;
entry.sigma_s2_used = sigma_s2;
entry.T_used = T;
end

function S = assign_result_entry(S,idx,entry,arrayName)
if isempty(S)
    if idx ~= 1
        error('The first assignment to %s must use index 1.',arrayName);
    end
    S = entry;
    return;
end
existingFields = sort(fieldnames(S));
newFields = sort(fieldnames(entry));
if ~isequal(existingFields,newFields)
    error('Field mismatch while assigning %s(%d).',arrayName,idx);
end
entry = orderfields(entry,S(1));
S(idx) = entry;
end

function T = make_empty_transition_diagnostic_table()
variableNames = { ...
    'Layer','CoarseIndex','PrevTargetDB','TargetDB','StepDB', ...
    'TargetKind','RefinementDepth','IsForced','IsInserted', ...
    'DualStateWindow','HandoffSource','BridgeAttempted', ...
    'BridgeNumAttempts','BridgeEtaUsed','BridgeTargetDB','BridgeCRBDB', ...
    'BridgeMinSINRDB','BridgeRelEigA','BridgeRelEigB','BridgeStatus', ...
    'BridgeAttemptSummary','ReportSuccess','ReportCRBDB', ...
    'ReportMinSINRDB','ReportRelEigA','ReportRelEigB','DeltaCRBDB', ...
    'Severity','NeedsProtection','RefinementInserted','InsertedMidDB', ...
    'BranchSwitchTried','RecoveryReason','SolutionSource', ...
    'RecoveryCovarianceGapRel','RecoveryStatus','MaxOuter','MaxPhaseInner', ...
    'RootCause','CarrierTargetCurrentMargin','CarrierTargetNextMargin', ...
    'CarrierActualNextMargin','CarrierCRBDeltaDB','CarrierQuality', ...
    'CarrierUsedPrevious','CarrierGoalFallback','CarrierBuildCVXStatus', ...
    'CarrierSelectedMode','CarrierOuterSteps','CarrierPhaseIOuterSteps', ...
    'CarrierAllPhaseSteps','CarrierHiddenAcceptedSteps', ...
    'CarrierGoalAttempts','CarrierBuildSeconds','ContinuationPath', ...
    'CarrierPath','CarrierPathSkipReason','AdditionalCandidatePaths', ...
    'TotalSolveSeconds'};
variableTypes = { ...
    'double','double','double','double','double', ...
    'cell','double','logical','logical', ...
    'logical','cell','logical', ...
    'double','double','double','double', ...
    'double','double','double','cell', ...
    'cell','logical','double', ...
    'double','double','double','double', ...
    'cell','logical','logical','double', ...
    'logical','cell','cell', ...
    'double','cell','double','double', ...
    'cell','double','double', ...
    'double','double','cell', ...
    'logical','logical','cell', ...
    'cell','double','double', ...
    'double','double', ...
    'double','double','cell', ...
    'cell','cell','cell', ...
    'double'};
T = table('Size',[0,numel(variableNames)], ...
    'VariableTypes',variableTypes,'VariableNames',variableNames);
end

function diag = make_default_transition_diag(gammaPrev,gammaNext,statusText)
if nargin < 3 || isempty(statusText)
    statusText = 'report-seed';
end
diag = struct( ...
    'attempted',false, ...
    'usedBridge',false, ...
    'handoffSource','report', ...
    'etaUsed',0, ...
    'bridgeTarget',gammaPrev, ...
    'bridgeTargetDB',10*log10(max(gammaPrev,realmin)), ...
    'bridgeCRBPhysicalDB',NaN, ...
    'minSINR',NaN, ...
    'minSINRdB',NaN, ...
    'relMinEigA',NaN, ...
    'relMinEigB',NaN, ...
    'status',char(statusText), ...
    'numAttempts',0, ...
    'gammaPrev',gammaPrev, ...
    'gammaNext',gammaNext, ...
    'attemptEtaList',zeros(1,0), ...
    'attemptTargetDB',zeros(1,0), ...
    'attemptAccepted',false(1,0), ...
    'attemptCRBPhysicalDB',zeros(1,0), ...
    'attemptMinSINRdB',zeros(1,0), ...
    'attemptRelMinEigA',zeros(1,0), ...
    'attemptRelMinEigB',zeros(1,0), ...
    'attemptStatus',{cell(1,0)}, ...
    'carrierInfo',struct());
end

function T = append_transition_diagnostic_row( ...
    T,L,coarseIndex,previousRecord,bestCandidate,task,useDualStateHere, ...
    transitionDiag,severity,transition,refinementInserted,insertedMid, ...
    branchSwitchTried,record,recoveryReason,Tused)

prevDB = NaN;
if isstruct(previousRecord) && isfield(previousRecord,'targetGammaDB')
    prevDB = previousRecord.targetGammaDB;
end
targetDB = 10*log10(max(task.gamma,realmin));
stepDB = targetDB-prevDB;
if task.isInsertedPoint
    targetKind = 'adaptive-refine';
elseif task.isForcedPoint
    targetKind = 'forced-grid';
else
    targetKind = 'formal';
end

reportSuccess = false;
reportCRBDB = NaN;
reportMinSINRDB = NaN;
reportRelA = NaN;
reportRelB = NaN;
if ~isempty(bestCandidate) && isstruct(bestCandidate) && ...
        isfield(bestCandidate,'sol') && is_solution_success(bestCandidate.sol)
    reportSuccess = true;
    solTmp = bestCandidate.sol;
    reportCRBDB = 10*log10(max( ...
        (solTmp.sigma_s2/max(Tused,1))*solTmp.metrics.CRB,realmin));
    reportMinSINRDB = 10*log10(max(min(solTmp.metrics.sinr),realmin));
    reportRelA = solTmp.metrics.crbInfo.relMinEigA;
    reportRelB = solTmp.metrics.crbInfo.relMinEigB;
end

solutionSource = 'pending-refinement';
recoveryGap = NaN;
recoveryStatus = 'not-retained';
if ~isempty(record) && isstruct(record)
    if isfield(record,'success') && logical(record.success)
        reportSuccess = true;
        reportCRBDB = record.CRBPhysicalDB;
        reportMinSINRDB = record.minSINRdB;
        reportRelA = record.relMinEigA;
        reportRelB = record.relMinEigB;
    end
    if isfield(record,'solutionSource') && ~isempty(record.solutionSource)
        solutionSource = record.solutionSource;
    end
    if isfield(record,'recoveryCovarianceGapRel')
        recoveryGap = record.recoveryCovarianceGapRel;
    end
    if isfield(record,'recoveryStatus') && ~isempty(record.recoveryStatus)
        recoveryStatus = record.recoveryStatus;
    end
end

attemptSummary = transition_attempt_summary(transitionDiag);
bridgeStatus = char(transitionDiag.status);
handoffSource = char(transitionDiag.handoffSource);
needsProtection = ~strcmp(char(severity),'none');
insertedMidDB = NaN;
if refinementInserted && isfinite(insertedMid)
    insertedMidDB = 10*log10(max(insertedMid,realmin));
end

maxOuter = NaN;
maxPhaseInner = NaN;
if ~isempty(record) && isstruct(record)
    if isfield(record,'maxOuterUsed'), maxOuter = record.maxOuterUsed; end
    if isfield(record,'maxPhaseInnerUsed'), maxPhaseInner = record.maxPhaseInnerUsed; end
end

carrierInfo=transitionDiag.carrierInfo;
buildInfo=struct();
if isstruct(carrierInfo)&&isfield(carrierInfo,'buildDetails')&& ...
        isstruct(carrierInfo.buildDetails)
    buildInfo=carrierInfo.buildDetails;
end
primaryPath=empty_ao_path_diagnostics('continuation');
carrierPath=empty_ao_path_diagnostics('commCarrier');
if isstruct(carrierInfo)&&isfield(carrierInfo,'primaryPath')
    primaryPath=carrierInfo.primaryPath;
end
if isstruct(carrierInfo)&&isfield(carrierInfo,'carrierPath')
    carrierPath=carrierInfo.carrierPath;
end
additionalPaths=[];
if isstruct(carrierInfo)&&isfield(carrierInfo,'additionalCandidates')
    additionalPaths=carrierInfo.additionalCandidates;
end
rootCause=classify_transition_root_cause( ...
    severity,transition,carrierInfo,primaryPath,carrierPath,record);
totalSolveSeconds=scalar_field(carrierInfo,'buildElapsedSeconds',0)+ ...
    scalar_field(primaryPath,'elapsedSeconds',0)+ ...
    scalar_field(carrierPath,'elapsedSeconds',0);

row = { ...
    L,coarseIndex,prevDB,targetDB,stepDB, ...
    {targetKind},task.depth,logical(task.isForcedPoint), ...
    logical(task.isInsertedPoint),logical(useDualStateHere), ...
    {handoffSource},logical(transitionDiag.attempted), ...
    transitionDiag.numAttempts,transitionDiag.etaUsed, ...
    transitionDiag.bridgeTargetDB,transitionDiag.bridgeCRBPhysicalDB, ...
    transitionDiag.minSINRdB,transitionDiag.relMinEigA, ...
    transitionDiag.relMinEigB,{bridgeStatus},{attemptSummary}, ...
    logical(reportSuccess),reportCRBDB,reportMinSINRDB,reportRelA,reportRelB, ...
    transition.deltaCRBDB,{char(severity)},logical(needsProtection), ...
    logical(refinementInserted),insertedMidDB,logical(branchSwitchTried), ...
    {char(recoveryReason)},{solutionSource},recoveryGap,{recoveryStatus}, ...
    maxOuter,maxPhaseInner, ...
    {rootCause},scalar_field(carrierInfo,'targetCurrentMargin',NaN), ...
    scalar_field(carrierInfo,'targetNextMargin',NaN), ...
    scalar_field(carrierInfo,'nextMargin',NaN), ...
    scalar_field(carrierInfo,'crbDeltaDB',NaN), ...
    {status_text(carrierInfo,'quality','none')}, ...
    logical_scalar_field(carrierInfo,'usedPreviousCarrier',false), ...
    logical_scalar_field(buildInfo,'goalFallbackUsed',false), ...
    {status_text(buildInfo,'lastCVXStatus','not-run')}, ...
    {status_text(buildInfo,'selectedMode','none')}, ...
    scalar_field(buildInfo,'outerAcceptedSteps',0), ...
    scalar_field(buildInfo,'phaseIOuterAcceptedSteps',0), ...
    scalar_field(buildInfo,'allPhaseAcceptedSteps',0), ...
    scalar_field(buildInfo,'hiddenAcceptedSteps',0), ...
    scalar_field(buildInfo,'goalAttemptCount',0), ...
    scalar_field(carrierInfo,'buildElapsedSeconds',0), ...
    {primaryPath},{carrierPath}, ...
    {status_text(carrierInfo,'carrierPathSkipReason','not-applicable')}, ...
    {additionalPaths},totalSolveSeconds};
T(end+1,:) = row;
end

function cause=classify_transition_root_cause( ...
    severity,transition,carrierInfo,primaryPath,carrierPath,record)
if strcmp(char(severity),'none')&& ...
        ~(isfinite(transition.deltaCRBDB)&&transition.deltaCRBDB<-0.10)
    if logical_scalar_field(carrierInfo,'selected',false)
        cause='carrier-improved-bestCRB';
    else
        cause='normal-continuation';
    end
    return;
end
primaryFirstBad=primaryPath.attempted&&~primaryPath.firstSDRFeasible;
carrierFirstBad=carrierPath.attempted&&~carrierPath.firstSDRFeasible;
if strcmp(char(severity),'failure')
    if primaryFirstBad&&carrierFirstBad
        cause='both-path-first-SDR-infeasible';
    elseif primaryFirstBad&&~carrierPath.attempted
        cause='first-SDR-infeasible-carrier-unavailable';
    elseif ~primaryPath.success&&carrierPath.attempted&&~carrierPath.success
        cause='both-AO-paths-failed';
    else
        cause='no-complete-feasible-candidate';
    end
    return;
end
if isfinite(transition.deltaCRBDB)&&transition.deltaCRBDB<-0.10
    cause='late-lower-CRB-branch-discovered';
    return;
end
selectedPath=primaryPath;
if ~isempty(record)&&isstruct(record)&& ...
        strcmp(status_text(record,'selectedInitialization',''), 'commCarrier')
    selectedPath=carrierPath;
elseif logical_scalar_field(carrierInfo,'selected',false)
    selectedPath=carrierPath;
end
quality=status_text(carrierInfo,'quality','none');
if strcmp(quality,'fallback')&& ...
        scalar_field(carrierInfo,'crbDeltaDB',0)>15
    cause='carrier-feasible-only-with-large-CRB-loss';
elseif ~carrierPath.attempted&&logical_scalar_field( ...
        carrierInfo,'attempted',false)
    cause='carrier-unavailable-or-identical';
elseif primaryFirstBad&&carrierPath.success
    cause='bestCRB-first-SDR-infeasible-carrier-rescue';
elseif selectedPath.phaseAccepted==0&&selectedPath.phaseFeasibilityRejects>0
    cause='inner-phase-blocked-by-SINR-feasibility';
elseif selectedPath.phaseAccepted==0&&selectedPath.phaseConditionRejects>0
    cause='inner-phase-blocked-by-FIM-conditioning';
elseif selectedPath.phaseAccepted==0&&selectedPath.phaseObjectiveRejects>0
    cause='inner-phase-blocked-by-CRB-acceptance';
elseif selectedPath.cvxFailureCount>0
    cause='CVX-failures-or-infeasible-fixed-phase';
elseif selectedPath.hitOuterBudget&&~selectedPath.outerStopTriggered
    cause='AO-iteration-budget-exhausted-before-convergence';
elseif primaryPath.success&&carrierPath.success
    cause='both-feasible-paths-converged-to-worse-CRB';
else
    cause='nonconvex-branch-transition';
end
end

function text = transition_attempt_summary(diag)
if ~diag.attempted || diag.numAttempts <= 0
    text = 'no bridge attempt';
    return;
end
recordedAttempts=min([numel(diag.attemptEtaList), ...
    numel(diag.attemptTargetDB),numel(diag.attemptCRBPhysicalDB), ...
    numel(diag.attemptMinSINRdB),numel(diag.attemptStatus)]);
if recordedAttempts<=0
    text=sprintf('%d Phase-I attempts; no per-attempt snapshot', ...
        diag.numAttempts);
    return;
end
parts = cell(1,recordedAttempts);
for ii = 1:recordedAttempts
    eta = diag.attemptEtaList(ii);
    targetDB = diag.attemptTargetDB(ii);
    crbDB = diag.attemptCRBPhysicalDB(ii);
    sinrDB = diag.attemptMinSINRdB(ii);
    status = diag.attemptStatus{ii};
    parts{ii} = sprintf('eta=%.2f,target=%.3f,CRB=%.3f,SINR=%.3f,%s', ...
        eta,targetDB,crbDB,sinrDB,status);
end
text = strjoin(parts,' | ');
if diag.numAttempts>recordedAttempts
    text=sprintf('%s [%d total Phase-I attempts]',text,diag.numAttempts);
end
end

function save_transition_diagnostics(fileName,T,scriptVersion,dualState, ...
    refine,diagnosticOptions)
if ~diagnosticOptions.enable
    return;
end
transitionDiagnostics = T; %#ok<NASGU>
save(fileName,'transitionDiagnostics','scriptVersion','dualState', ...
    'refine','diagnosticOptions','-v7.3');
append_compact_transition_trace(T,diagnosticOptions);
end

function initialize_compact_transition_trace(options,scriptVersion)
if ~logical_scalar_field(options,'writeCompactPointTrace',false)
    return;
end
fileName=status_text(options,'traceTextFile','');
if isempty(fileName),return;end
fid=fopen(fileName,'w');
if fid<0
    warning('SIM_CRB_SINR:PointTraceOpenFailed', ...
        'Cannot initialize compact point diagnostics: %s',fileName);
    return;
end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'SIM persistent bestCRB/commCarrier root-cause trace\n');
fprintf(fid,'Version: %s\n',scriptVersion);
fprintf(fid,['Each POINT line records both AO paths; DETAIL lines are ', ...
    'emitted only for abnormal transitions.\n\n']);
end

function append_compact_transition_trace(T,options)
persistent lastFile lastRowCount
if isempty(T)||~logical_scalar_field(options,'writeCompactPointTrace',false)
    return;
end
fileName=status_text(options,'traceTextFile','');
if isempty(fileName),return;end
if isempty(lastFile)||~strcmp(lastFile,fileName)|| ...
        isempty(lastRowCount)||height(T)<lastRowCount
    lastFile=fileName;
    lastRowCount=0;
end
if height(T)<=lastRowCount,return;end
fid=fopen(fileName,'a');
if fid<0
    warning('SIM_CRB_SINR:PointTraceAppendFailed', ...
        'Cannot append compact point diagnostics: %s',fileName);
    return;
end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
for ii=lastRowCount+1:height(T)
    primary=T.ContinuationPath{ii};
    carrier=T.CarrierPath{ii};
    fprintf(fid,['POINT L=%d target=%.5f prev=%.5f depth=%d ', ...
        'dCRB=%+.3f severity=%s source=%s root=%s ', ...
        'reserve=%+.3e/%.3e quality=%s inherited=%d ', ...
        'best=%s:%.3f:%s carrier=%s:%.3f:%s time=%.1fs\n'], ...
        T.Layer(ii),T.TargetDB(ii),T.PrevTargetDB(ii), ...
        T.RefinementDepth(ii),T.DeltaCRBDB(ii), ...
        T.Severity{ii},T.SolutionSource{ii},T.RootCause{ii}, ...
        T.CarrierActualNextMargin(ii),T.CarrierTargetNextMargin(ii), ...
        T.CarrierQuality{ii},T.CarrierUsedPrevious(ii), ...
        compact_path_status(primary),primary.physicalCRBDB, ...
        primary.firstCVXStatus,compact_path_status(carrier), ...
        carrier.physicalCRBDB,carrier.firstCVXStatus, ...
        T.TotalSolveSeconds(ii));
    abnormal=~strcmp(T.Severity{ii},'none')|| ...
        T.RefinementInserted(ii)||T.BranchSwitchTried(ii)|| ...
        strcmp(T.CarrierQuality{ii},'fallback')|| ...
        (isfinite(T.DeltaCRBDB(ii))&&T.DeltaCRBDB(ii)<-0.10);
    if abnormal
        fprintf(fid,['DETAIL firstOK=%d/%d outer=%d/%d ', ...
            'Wacc=%d/%d PhiAcc=%d/%d inner=%d/%d ', ...
            'reject(feas,obj,cond,model)=%d,%d,%d,%d/', ...
            '%d,%d,%d,%d CVXfail=%d/%d retry=%d/%d ', ...
            'carrierBuildCVX=%s carrierSteps(out,phaseI,all)=%d,%d,%d ', ...
            'goalFallback=%d carrierSkip=%s reason=%s\n'], ...
            primary.firstSDRFeasible,carrier.firstSDRFeasible, ...
            primary.outerIterations,carrier.outerIterations, ...
            primary.wAccepted,carrier.wAccepted, ...
            primary.phaseAccepted,carrier.phaseAccepted, ...
            primary.innerPhaseSteps,carrier.innerPhaseSteps, ...
            primary.phaseFeasibilityRejects, ...
            primary.phaseObjectiveRejects,primary.phaseConditionRejects, ...
            primary.phaseModelRejects,carrier.phaseFeasibilityRejects, ...
            carrier.phaseObjectiveRejects,carrier.phaseConditionRejects, ...
            carrier.phaseModelRejects,primary.cvxFailureCount, ...
            carrier.cvxFailureCount,primary.sdrRetryCount, ...
            carrier.sdrRetryCount,T.CarrierBuildCVXStatus{ii}, ...
            T.CarrierOuterSteps(ii),T.CarrierPhaseIOuterSteps(ii), ...
            T.CarrierAllPhaseSteps(ii),T.CarrierGoalFallback(ii), ...
            T.CarrierPathSkipReason{ii},T.RecoveryReason{ii});
        if primary.hitOuterBudget||carrier.hitOuterBudget
            fprintf(fid,['BUDGET outer=%d/%g,%d/%g stop=%d/%d ', ...
                'lastDecrease=%.3e/%.3e grad=%.3e/%.3e\n'], ...
                primary.outerIterations,primary.outerIterationBudget, ...
                carrier.outerIterations,carrier.outerIterationBudget, ...
                primary.outerStopTriggered,carrier.outerStopTriggered, ...
                primary.lastRelativeDecrease,carrier.lastRelativeDecrease, ...
                primary.lastPhaseGradNorm,carrier.lastPhaseGradNorm);
        end
    end
end
lastRowCount=height(T);
end

function run_transition_diagnostic_selftest(packageDir,scriptVersion)
resultsDir=fullfile(packageDir,'results');
if exist(resultsDir,'dir')~=7,mkdir(resultsDir);end
options=struct('enable',true,'writeCompactPointTrace',true, ...
    'traceTextFile',fullfile(resultsDir, ...
    'SIM_CRB_SINR_point_diagnostics_selftest.txt'));
diagnosticsFile=fullfile(resultsDir, ...
    'SIM_CRB_SINR_diagnostics_selftest.mat');
initialize_compact_transition_trace(options,scriptVersion);
T=make_empty_transition_diagnostic_table();
previous=struct('targetGammaDB',10);
task=struct('gamma',11,'depth',2,'isInsertedPoint',false, ...
    'isForcedPoint',true);
transition=struct('deltaCRBDB',12.5);
diag=make_default_transition_diag(10,11,'selftest-carrier');
diag.attempted=true;
diag.handoffSource='commCarrier';
diag.numAttempts=4;
diag.attemptEtaList=0.25;
diag.attemptTargetDB=10*log10(11);
diag.attemptCRBPhysicalDB=-20;
diag.attemptMinSINRdB=10*log10(11);
diag.attemptStatus={'Solved'};
failureHist=struct('failureStage','first-SDR', ...
    'firstSDRFeasible',false,'firstSDRStatus','Infeasible', ...
    'finalSDRStatus','Infeasible','outerIterationsCompleted',1, ...
    'initialIncumbentFeasibleAtTarget',false, ...
    'initialIncumbentMinSINR',9,'sdrRetryCount',2, ...
    'finalSDRRetryCount',1,'WBlockAccepted',false, ...
    'phaseBlockAccepted',false,'phaseInner',0, ...
    'outerMarginAcceptedSteps',1,'phaseFeasibilityRejects',3, ...
    'phaseObjectiveRejects',0,'phaseConditionRejects',0, ...
    'phaseModelRejects',1,'phaseLineSearchTrials',4, ...
    'phaseGradNorm',0.5,'sdrCVXStatus',{{'Infeasible'}});
failureSol=struct('success',false,'status','CVX infeasible');
primary=summarize_ao_path('continuation',failureSol,failureHist,11,10,2);
successHist=failureHist;
successHist.failureStage='none';
successHist.firstSDRFeasible=true;
successHist.firstSDRStatus='Solved';
successHist.finalSDRStatus='Solved';
successHist.phaseBlockAccepted=true;
successHist.phaseFeasibilityRejects=0;
successHist.sdrCVXStatus={'Solved','Inaccurate/Solved'};
successSol=struct('success',true,'status','Solved', ...
    'sigma_s2',1,'bestIteration',1,'bestStage','phase', ...
    'metrics',struct('CRB',0.1,'crbValid',true,'sinrFeasible',true, ...
    'sinr',[11.2;11.4], ...
    'crbInfo',struct('relMinEigA',1e-4,'relMinEigB',2e-4)));
carrier=summarize_ao_path('commCarrier',successSol,successHist,11,10,3);
assert(primary.cvxFailureCount==1&&carrier.cvxInaccurateCount==1&& ...
    carrier.phaseAccepted==1,'AO history summaries lost CVX/phase details.');
diag.carrierInfo=struct('attempted',true,'selected',true, ...
    'targetCurrentMargin',0.10,'targetNextMargin',0.05, ...
    'nextMargin',0.02,'crbDeltaDB',18,'quality','fallback', ...
    'usedPreviousCarrier',true,'buildElapsedSeconds',1, ...
    'primaryPath',primary,'carrierPath',carrier, ...
    'carrierPathSkipReason','','additionalCandidates',[], ...
    'buildDetails',struct('goalFallbackUsed',true, ...
    'lastCVXStatus','Solved','selectedMode','phase-I', ...
    'outerAcceptedSteps',2,'phaseIOuterAcceptedSteps',1, ...
    'allPhaseAcceptedSteps',1,'hiddenAcceptedSteps',1, ...
    'goalAttemptCount',2));
record=struct('success',true,'CRBPhysicalDB',-20, ...
    'minSINRdB',10*log10(11),'relMinEigA',1e-4,'relMinEigB',2e-4, ...
    'solutionSource','persistent-commCarrier', ...
    'selectedInitialization','commCarrier', ...
    'maxOuterUsed',20,'maxPhaseInnerUsed',4);
T=append_transition_diagnostic_row(T,3,7,previous,[],task,true, ...
    diag,'hard',transition,false,NaN,false,record,'selftest-jump',10);
assert(height(T)==1&&width(T)>50,'Diagnostic table schema is incomplete.');
assert(contains(T.BridgeAttemptSummary{1},'4 total Phase-I attempts'), ...
    'Phase-I attempt summary did not safely capture all attempts.');
assert(strcmp(T.RootCause{1},'carrier-feasible-only-with-large-CRB-loss'), ...
    'Jump root-cause classification failed.');
save_transition_diagnostics(diagnosticsFile,T,scriptVersion, ...
    struct(),struct(),options);
saved=load(diagnosticsFile,'transitionDiagnostics');
assert(height(saved.transitionDiagnostics)==1, ...
    'Diagnostic MAT did not retain the target record.');
trace=fileread(options.traceTextFile);
assert(contains(trace,'POINT L=3')&&contains(trace,'DETAIL firstOK=0/1'), ...
    'Compact TXT did not retain both point and abnormal details.');
reportFile=fullfile(resultsDir,'SIM_CRB_SINR_diagnostics_selftest_report.txt');
SIM_print_transition_diagnostics('summary',diagnosticsFile, ...
    fullfile(resultsDir,'selftest-no-formal-results.mat'),reportFile);
summaryText=fileread(reportFile);
assert(contains(summaryText,'bestCRB AO attempted / solved')&& ...
    contains(summaryText,'carrier-feasible-only-with-large-CRB-loss'), ...
    'Compact final summary omitted dual-path root-cause evidence.');
fprintf(['DIAGNOSTIC_SELFTEST_OK columns=%d rows=%d root=%s ', ...
    'traceBytes=%d reportBytes=%d\nMAT=%s\nTXT=%s\nREPORT=%s\n'], ...
    width(T),height(T),T.RootCause{1},numel(trace),numel(summaryText), ...
    diagnosticsFile,options.traceTextFile,reportFile);
end

function value = scalar_field(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && isscalar(s.(name)) ...
        && ~isempty(s.(name))
    value = s.(name);
end
end

function value = logical_scalar_field(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && isscalar(s.(name)) ...
        && ~isempty(s.(name))
    value = logical(s.(name));
end
end

function textValue = status_text(s,fieldName,defaultValue)
if nargin < 3
    defaultValue = 'not-reported';
end
textValue = defaultValue;
if isstruct(s) && isfield(s,fieldName) && ~isempty(s.(fieldName))
    textValue = char(string(s.(fieldName)));
end
end
