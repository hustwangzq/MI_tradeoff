function params = default_params()
%DEFAULT_PARAMS Default parameters for the MI-MMSE RIS/SIM ISAC demo.
% This file only sets simulation parameters. There is no direct BS-to-user
% or BS-to-target link in the code. All communication and sensing paths pass
% through the programmable metasurface architecture.

params.c0 = 3e8;                                  % Speed of light.
params.fc = 5.8e9;                                % Carrier frequency.
params.lambda = params.c0/params.fc;              % Wavelength.
params.d = params.lambda/2;                       % Element spacing.

% --------------------------- Array dimensions ---------------------------
params.NtSide = 3;                                % BS UPA side length; Nt=NtSide^2.
params.Nt = params.NtSide^2;                      % Number of BS antennas.
params.Nr = 3;                                    % RIS/SIM rows used in this demo.
params.Nc = 3;                                    % RIS/SIM columns used in this demo.
params.N = params.Nr*params.Nc;                   % Elements per layer, N=Nr*Nc.
params.K = 2;                                     % Number of communication users.
params.Qtar = 2;                                  % Number of point targets.
params.T = 16;                                    % Number of sensing snapshots.
params.Ns = params.K + params.Nt;                 % Total streams: K communication + Nt sensing.

% ---------------------- Power, weights, and noise ------------------------
params.P0 = 0.2;                                  % BS transmit power budget.
params.wc = 0.5;                                  % Default communication MI weight.
params.ws = 0.5;                                  % Default sensing MI weight.
params.userWeight = ones(params.K,1);             % Per-user weights; set to one to match the paper.
params.SNRc_dB = 15;                              % Reference communication SNR for sigma_c^2.
params.SNRs_dB = 15;                              % Reference sensing SNR for sigma_s^2.
params.noise.roundToOrder = true;                 % Report clean powers such as 1e-10.
params.noise.minSigmaC2 = 1e-16;
params.noise.minSigmaS2 = 1e-18;

% -------------------------- Tradeoff settings ---------------------------
params.tradeoff.rhoGrid = 0.001:0.09:0.991;          % wc=rho, ws=1-rho for MI tradeoff points.
params.tradeoff.rhoConv = 0.50;                   % Weight used only for convergence curves.

% --------------------------- Channel model ------------------------------
params.channel.model = 'Rician';                  % Far-field segments use Rician fading.
params.ricianK = 5;                               % Linear Rician factor kappa, not dB.
params.pathlossRef_dB = -30;                      % Pathloss at 1 m for far-field segments.
params.pathlossExp = 2.2;                         % Pathloss exponent for far-field segments.

% ---------------------- Target-response parameters ----------------------
params.target.sigmaAlpha2 = 1;                    % Reflection coefficient variance before pathloss scaling.
params.target.rhDiagLoad = 1e-8;                  % Tiny background loading for low-rank R_hs.

% -------------------------- SIM geometry --------------------------------
params.sim.totalThickness = 3*params.lambda;      % BS-to-outermost-SIM-layer distance = 3 lambda.
params.sim.userDistance = 2.4;                    % Actual outermost layer-to-user distance.
params.sim.targetDistance = 2.6;                  % Actual outermost layer-to-target distance.
params.sim.userAzDeg = linspace(-5,5,params.K); % User azimuth angles relative to outermost-layer broadside.
params.sim.userElDeg = zeros(1,params.K);         % User elevation angles.
params.sim.targetAzDeg = linspace(30,40,params.Qtar); % Target azimuths; separated from users.
params.sim.targetElDeg = linspace(-4,4,params.Qtar);  % Target elevations.

% ------------------------ Multi-hop rRIS geometry -----------------------
params.rris.hopDistance = 0.8;                    % Inter-rRIS distance; should exceed Rayleigh boundary.
params.rris.lastHopUserDistance = 0.9;            % Actual last-rRIS-to-user distance.
params.rris.lastHopTargetDistance = 0.8;          % Actual last-rRIS-to-target distance.
params.rris.userAzDeg = linspace(-5,5,params.K);% User azimuth angles after the last rRIS.
params.rris.userElDeg = zeros(1,params.K);        % User elevation angles.
params.rris.targetAzDeg = linspace(30,40,params.Qtar); % Target azimuths after the last rRIS.
params.rris.targetElDeg = linspace(-4,4,params.Qtar);  % Target elevations.

% -------------------------- Random seeds --------------------------------
params.seed.global = 7;                           % Global seed used by main scripts.
params.seed.init = 1000;                          % Initialization seed for W, theta, and S.
params.seed.userBase = 2100;                      % Seed base for user-link NLoS components.
params.seed.targetBase = 3100;                    % Seed base for target response coefficients.
params.seed.rrisBase = 4100;                      % Seed base for rRIS far-field hop NLoS components.

% ------------------------ Algorithm parameters --------------------------
params.alg.maxOuter = 800;                        % Maximum outer BCD iterations.
params.alg.maxPhaseInner = 100;                    % Maximum inner phase-update iterations.
params.alg.useOuterStop = false;                  % For convergence plots, run fixed iterations by default.
params.alg.minOuter = 30;                         % Minimum iterations before optional stopping.
params.alg.tolOuter = 0;                       % Optional relative convergence threshold for WSMI.
params.alg.tolPhase = 0;                       % Optional inner phase objective threshold.
params.alg.tolLambda = 1e-16;                      % Relative accuracy for lambda bisection.
params.alg.maxBisect = 100;                        % Maximum bisection steps for lambda.
params.alg.armijo_alpha0 = 1;                     % Initial Armijo step size.
params.alg.armijo_beta = 0.5;                     % Armijo step shrinking factor.
params.alg.armijo_c = 1e-4;                       % Armijo sufficient-decrease parameter.
params.alg.armijo_minAlpha = 1e-12;               % Minimum line-search step.
params.alg.phaseGradTol = 0;                   % Stop phase inner loop if gradient is tiny.
params.alg.minRcondGamma = 1e-18;                 % Numerical check for lambda=0 linear solve.
params.alg.hessianReg = 0;                        % Default zero; lambda, not regularization, handles power constraint.
params.alg.verbose = false;                       % Print per-iteration information if true.
params.alg.debugBlock = true;                     % Print debugblock information if true.
params.alg.verbose = true;

% Feasible monotone proximal-SCA phase-update defaults. The dedicated
% SIM_CRB_SINR.m script overrides the iteration caps with a balanced test
% profile, but these fields keep other callers compatible.
params.alg.phaseMethod = 'FP-SCA';
params.alg.outerMinIterations = 15;
params.alg.outerStableIterations = 6;
params.alg.phaseExactConstraintTol = 2e-5;
params.alg.phaseSubproblemConstraintTol = 2e-5;
params.alg.phaseTrustRadius0 = 0.08;
params.alg.phaseTrustRadiusMin = 1e-5;
params.alg.phaseTrustRadiusMax = 0.20;
params.alg.phaseTrustShrink = 0.5;
params.alg.phaseTrustExpand = 1.10;
params.alg.phaseCurvatureF0 = 1;
params.alg.phaseCurvatureC0 = 1;
params.alg.phaseCurvatureGrowth = 2;
params.alg.phaseCurvatureDecrease = 1.10;
params.alg.phaseCurvatureMin = 1e-4;
params.alg.phaseCurvatureMax = 1e8;
params.alg.phaseMaxModelResolve = 2;
params.alg.phaseMaxLineSearch = 8;
params.alg.phaseLineSearchBeta = 0.5;
params.alg.phaseMinStepScale = 1/128;
params.alg.phaseStepTol = 1e-6;
params.alg.phaseModelTol = 1e-7;
params.alg.phaseObjectiveTol = 1e-9;
params.alg.phaseSufficientDecrease = 1e-4;
params.alg.phaseRequireModelUpperBound = false;
end
