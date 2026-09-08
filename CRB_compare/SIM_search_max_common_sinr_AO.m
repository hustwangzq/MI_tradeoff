function [gammaMaxFound,solAtGammaMax,searchInfo] = ...
    SIM_search_max_common_sinr_AO(params,ch,alg,initState,gammaStart,searchOpts)
%SEARCH_MAX_COMMON_SINR_AO_SIM SIM-specific common-SINR range search.
%
% This function is intentionally separate from SIM_search_max_common_sinr_AO.m
% so that the already validated SIM workflow is not changed.
%
% Supported continuation modes:
%   multiplicative : Gamma_next = growthFactor * Gamma_success.
%                    Intended for the one-layer SIM search.
%   additive       : Gamma_next = Gamma_success + currentStep.
%                    Intended for the two-/three-layer SIM searches. When a
%                    trial fails, currentStep is reduced and the search is
%                    retried from the latest verified feasible point.
%
% Optional feasibility recovery:
%   For L=2/3, gammaStart may be inherited from the previous layer and may
%   not be immediately feasible because the physical propagation structure
%   changes with L. If enableRecovery=true, the function reduces gammaStart
%   by recoveryFactor until a feasible current-layer starting point is found.
%
% The returned gammaMaxFound is the largest SINR value verified by this AO
% continuation under the supplied initialization and iteration budgets. It
% is not claimed to be the global max-min SINR optimum.
%
% Optional full-budget boundary confirmation:
%   After the normal continuation/refinement, the current maximum is rerun
%   with alg.maxOuter. If confirmed, small additive boundary steps test
%   whether the short-budget range search stopped prematurely.

if nargin < 6 || isempty(searchOpts)
    searchOpts = struct();
end
searchOpts = fill_defaults(searchOpts);

validateattributes(gammaStart,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'gammaStart');
validateattributes(searchOpts.maxOuterGrowth,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});
validateattributes(searchOpts.maxOuterConfirm,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});
validateattributes(searchOpts.maxOuterRefine,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});
validateattributes(searchOpts.maxExpand,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});
validateattributes(searchOpts.maxRefine,{'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative'});
validateattributes(searchOpts.rangeRelTol,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(searchOpts.recoveryFactor,{'numeric'}, ...
    {'real','finite','scalar','>',0,'<',1});
validateattributes(searchOpts.maxRecovery,{'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative'});
validateattributes(searchOpts.enableFullBudgetBoundary, ...
    {'logical','numeric'},{'scalar'});
validateattributes(searchOpts.maxBoundaryTrials,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});
validateattributes(searchOpts.boundaryStepRel,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(searchOpts.boundaryStepMinRel,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(searchOpts.boundaryStepMinAbs,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(searchOpts.boundaryStepShrink,{'numeric'}, ...
    {'real','finite','scalar','>',0,'<',1});

searchMode = lower(char(string(searchOpts.searchMode)));
if ~ismember(searchMode,{'multiplicative','additive'})
    error('SIM_search_max_common_sinr_AO:InvalidSearchMode', ...
        'searchMode must be ''multiplicative'' or ''additive''.');
end

if strcmp(searchMode,'multiplicative')
    validateattributes(searchOpts.growthFactor,{'numeric'}, ...
        {'real','finite','scalar','>',1},mfilename,'growthFactor');
else
    validateattributes(searchOpts.additiveStep,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'additiveStep');
    validateattributes(searchOpts.additiveStepMin,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'additiveStepMin');
    validateattributes(searchOpts.stepShrink,{'numeric'}, ...
        {'real','finite','scalar','>',0,'<',1},mfilename,'stepShrink');
    validateattributes(searchOpts.maxStepHalving,{'numeric'}, ...
        {'real','finite','scalar','integer','positive'});
end

trace = struct([]);
trialCounter = 0;
lastSol = make_empty_failure_solution(gammaStart);
lastHist = struct();
if isfield(ch,'type') && strcmpi(ch.type,'rRIS')
    rangeTag = 'rRIS';
else
    rangeTag = 'SIM';
end

fprintf(['    [%s-range] requested Gamma=%.6e (%.3f dB), mode=%s, ', ...
    'budgets=%d/%d/%d\n'],rangeTag,gammaStart,10*log10(gammaStart),searchMode, ...
    searchOpts.maxOuterGrowth,searchOpts.maxOuterConfirm, ...
    searchOpts.maxOuterRefine);

% -------------------------------------------------------------------------
% 1) Find a feasible current-layer starting point.
% -------------------------------------------------------------------------
gammaTry = gammaStart;
startFound = false;
recoveryUsed = 0;

for ir = 0:searchOpts.maxRecovery
    if ir == 0
        stageGrowth = 'start-growth';
        stageConfirm = 'start-confirm';
    else
        stageGrowth = 'recovery-growth';
        stageConfirm = 'recovery-confirm';
    end

    fprintf('    [%s-range] start try %02d: Gamma=%.6e (%.3f dB)\n', ...
        rangeTag,ir,gammaTry,10*log10(gammaTry));

    [solTry,histTry] = run_trial(params,ch,alg,initState,gammaTry, ...
        searchOpts.maxOuterGrowth);
    [trace,trialCounter] = add_record(trace,trialCounter,stageGrowth, ...
        gammaTry,searchOpts.maxOuterGrowth,solTry,NaN);
    lastSol = solTry;
    lastHist = histTry;

    if ~solTry.success && ~is_first_sdr_failure(solTry)
        [solConfirm,histConfirm] = run_trial(params,ch,alg,initState, ...
            gammaTry,searchOpts.maxOuterConfirm);
        [trace,trialCounter] = add_record(trace,trialCounter,stageConfirm, ...
            gammaTry,searchOpts.maxOuterConfirm,solConfirm,NaN);
        lastSol = solConfirm;
        lastHist = histConfirm;
        if solConfirm.success
            solTry = solConfirm;
            histTry = histConfirm;
        end
    end

    if solTry.success
        startFound = true;
        solStart = solTry;
        histStart = histTry;
        gammaStartFeasible = gammaTry;
        recoveryUsed = ir;
        break;
    end

    if ~searchOpts.enableRecovery
        break;
    end

    gammaNext = max(searchOpts.recoveryMinGamma, ...
        searchOpts.recoveryFactor*gammaTry);
    if gammaNext >= gammaTry*(1-10*eps) || ...
            gammaTry <= searchOpts.recoveryMinGamma*(1+10*eps)
        break;
    end
    gammaTry = gammaNext;
end

if ~startFound
    gammaMaxFound = NaN;
    solAtGammaMax = lastSol;
    searchInfo = make_failure_search_info(searchOpts,trace,gammaStart, ...
        lastHist,'No feasible current-layer starting point was found.');
    fprintf(['    [%s-range] failed to find a feasible starting point ', ...
        'from requested Gamma %.6e.\n'],rangeTag,gammaStart);
    return;
end

gammaSuccess = gammaStartFeasible;
solSuccess = solStart;
histSuccess = histStart;
gammaFail = NaN;
upperBoundFound = false;
expandUsed = 0;
confirmUsed = 0;
bridgeUsed = 0;
stepReductionUsed = 0;

% -------------------------------------------------------------------------
% 2) Continue upward using the selected mode.
% -------------------------------------------------------------------------
if strcmp(searchMode,'multiplicative')
    for ie = 1:searchOpts.maxExpand
        expandUsed = ie;
        gammaUpper = searchOpts.growthFactor*gammaSuccess;
        warmUpper = make_warm_init(initState,solSuccess);

        fprintf(['    [%s-range] multiply %03d: Gamma=%.6e ', ...
            '(%.3f dB)\n'],rangeTag,ie,gammaUpper,10*log10(gammaUpper));

        [solUpper,histUpper] = run_trial(params,ch,alg,warmUpper, ...
            gammaUpper,searchOpts.maxOuterGrowth);
        [trace,trialCounter] = add_record(trace,trialCounter, ...
            'multiply-growth',gammaUpper,searchOpts.maxOuterGrowth, ...
            solUpper,NaN);

        if ~solUpper.success && ~is_first_sdr_failure(solUpper)
            confirmUsed = confirmUsed + 1;
            [solConfirm,histConfirm] = run_trial(params,ch,alg,warmUpper, ...
                gammaUpper,searchOpts.maxOuterConfirm);
            [trace,trialCounter] = add_record(trace,trialCounter, ...
                'multiply-confirm',gammaUpper,searchOpts.maxOuterConfirm, ...
                solConfirm,NaN);
            if solConfirm.success
                solUpper = solConfirm;
                histUpper = histConfirm;
            end
        end

        if solUpper.success
            gammaSuccess = gammaUpper;
            solSuccess = solUpper;
            histSuccess = histUpper;
            continue;
        end

        % Bridge the failed multiplicative jump by geometric midpoints. A
        % successful bridge is used to retry the original failed upper point.
        gammaUpperCurrent = gammaUpper;
        upperRecovered = false;
        for ib = 1:searchOpts.maxBridge
            if gammaUpperCurrent/gammaSuccess - 1 <= searchOpts.rangeRelTol
                break;
            end

            bridgeUsed = bridgeUsed + 1;
            gammaBridge = sqrt(gammaSuccess*gammaUpperCurrent);
            warmBridge = make_warm_init(initState,solSuccess);

            fprintf(['    [%s-range] bridge %02d: Gamma=%.6e ', ...
                'between %.6e and %.6e\n'],rangeTag,ib,gammaBridge, ...
                gammaSuccess,gammaUpperCurrent);

            [solBridge,histBridge] = run_trial(params,ch,alg,warmBridge, ...
                gammaBridge,searchOpts.maxOuterRefine);
            [trace,trialCounter] = add_record(trace,trialCounter,'bridge', ...
                gammaBridge,searchOpts.maxOuterRefine,solBridge,NaN);

            if ~solBridge.success
                gammaUpperCurrent = gammaBridge;
                continue;
            end

            gammaSuccess = gammaBridge;
            solSuccess = solBridge;
            histSuccess = histBridge;

            warmRetry = make_warm_init(initState,solSuccess);
            [solRetry,histRetry] = run_trial(params,ch,alg,warmRetry, ...
                gammaUpperCurrent,searchOpts.maxOuterConfirm);
            [trace,trialCounter] = add_record(trace,trialCounter, ...
                'bridge-retry-upper',gammaUpperCurrent, ...
                searchOpts.maxOuterConfirm,solRetry,NaN);

            if solRetry.success
                gammaSuccess = gammaUpperCurrent;
                solSuccess = solRetry;
                histSuccess = histRetry;
                upperRecovered = true;
                break;
            end
        end

        if upperRecovered
            continue;
        end

        gammaFail = gammaUpperCurrent;
        upperBoundFound = true;
        break;
    end
else
    currentStep = searchOpts.additiveStep;
    stepHalvingStreak = 0;

    fprintf(['    [%s-range] additive step=%.6e, minimum step=%.6e, ', ...
        'shrink=%.3f\n'],rangeTag,currentStep,searchOpts.additiveStepMin, ...
        searchOpts.stepShrink);

    for ie = 1:searchOpts.maxExpand
        expandUsed = ie;
        gammaUpper = gammaSuccess + currentStep;
        warmUpper = make_warm_init(initState,solSuccess);

        fprintf(['    [%s-range] add %03d: Gamma=%.6e (%.3f dB), ', ...
            'step=%.6e\n'],rangeTag,ie,gammaUpper,10*log10(gammaUpper), ...
            currentStep);

        [solUpper,histUpper] = run_trial(params,ch,alg,warmUpper, ...
            gammaUpper,searchOpts.maxOuterGrowth);
        [trace,trialCounter] = add_record(trace,trialCounter, ...
            'add-growth',gammaUpper,searchOpts.maxOuterGrowth, ...
            solUpper,currentStep);

        if ~solUpper.success && ~is_first_sdr_failure(solUpper)
            confirmUsed = confirmUsed + 1;
            [solConfirm,histConfirm] = run_trial(params,ch,alg,warmUpper, ...
                gammaUpper,searchOpts.maxOuterConfirm);
            [trace,trialCounter] = add_record(trace,trialCounter, ...
                'add-confirm',gammaUpper,searchOpts.maxOuterConfirm, ...
                solConfirm,currentStep);
            if solConfirm.success
                solUpper = solConfirm;
                histUpper = histConfirm;
            end
        end

        if solUpper.success
            gammaSuccess = gammaUpper;
            solSuccess = solUpper;
            histSuccess = histUpper;
            stepHalvingStreak = 0;
            continue;
        end

        % The failed value is not treated as a permanent boundary yet.
        % Reduce the fixed linear increment and retry from gammaSuccess.
        gammaFail = gammaUpper;
        nextStep = currentStep*searchOpts.stepShrink;
        stepHalvingStreak = stepHalvingStreak + 1;
        stepReductionUsed = stepReductionUsed + 1;

        if nextStep < searchOpts.additiveStepMin || ...
                stepHalvingStreak >= searchOpts.maxStepHalving
            upperBoundFound = true;
            break;
        end

        currentStep = nextStep;
        fprintf(['    [%s-range] failed upper trial; reducing additive ', ...
            'step to %.6e and retrying from Gamma=%.6e.\n'], ...
            rangeTag,currentStep,gammaSuccess);
    end
end

if ~upperBoundFound
    warning('SIM_search_max_common_sinr_AO:NoFailedUpperBound', ...
        ['No failed upper trial was retained within maxExpand=%d. The ', ...
         'returned Gamma is the largest verified point, not a bracketed ', ...
         'boundary.'],searchOpts.maxExpand);
end

% -------------------------------------------------------------------------
% 3) Refine the final feasible/failed bracket.
% -------------------------------------------------------------------------
refineUsed = 0;
if upperBoundFound && isfinite(gammaFail) && gammaFail > gammaSuccess
    for ir = 1:searchOpts.maxRefine
        relativeGap = (gammaFail-gammaSuccess)/max(gammaSuccess,realmin);
        absoluteGap = gammaFail-gammaSuccess;

        if relativeGap <= searchOpts.rangeRelTol
            break;
        end
        if strcmp(searchMode,'additive') && ...
                absoluteGap <= searchOpts.additiveStepMin
            break;
        end

        if strcmp(searchMode,'multiplicative')
            gammaMid = sqrt(gammaSuccess*gammaFail);
        else
            gammaMid = 0.5*(gammaSuccess+gammaFail);
        end

        warmMid = make_warm_init(initState,solSuccess);
        fprintf(['    [%s-range] refine %02d: Gamma=%.6e, ', ...
            'relative gap=%.4f%%\n'],rangeTag,ir,gammaMid,100*relativeGap);

        [solMid,histMid] = run_trial(params,ch,alg,warmMid,gammaMid, ...
            searchOpts.maxOuterRefine);
        [trace,trialCounter] = add_record(trace,trialCounter,'refine', ...
            gammaMid,searchOpts.maxOuterRefine,solMid, ...
            gammaMid-gammaSuccess);
        refineUsed = ir;

        if solMid.success
            gammaSuccess = gammaMid;
            solSuccess = solMid;
            histSuccess = histMid;
        else
            gammaFail = gammaMid;
        end
    end
end

% -------------------------------------------------------------------------
% 4) Confirm and locally extend the boundary using the full AO budget.
% -------------------------------------------------------------------------
gammaBeforeBoundaryConfirm = gammaSuccess;
boundaryTrialsUsed = 0;
boundaryFailedUpper = NaN;
boundaryFinalStep = NaN;
fullBudgetConfirmed = false;

if logical(searchOpts.enableFullBudgetBoundary)
    fullBudget = alg.maxOuter;
    warmFull = make_warm_init(initState,solSuccess);

    fprintf(['    [%s-range] full-budget confirmation at Gamma=%.6e ', ...
        '(%.3f dB), maxOuter=%d\n'],rangeTag,gammaSuccess, ...
        10*log10(gammaSuccess),fullBudget);

    [solFull,histFull] = run_trial(params,ch,alg,warmFull, ...
        gammaSuccess,fullBudget);
    boundaryTrialsUsed = boundaryTrialsUsed + 1;
    [trace,trialCounter] = add_record(trace,trialCounter, ...
        'boundary-confirm-current',gammaSuccess,fullBudget,solFull,0);

    if solFull.success
        fullBudgetConfirmed = true;
        solSuccess = solFull;
        histSuccess = histFull;

        boundaryStepMin = max(searchOpts.boundaryStepMinAbs, ...
            searchOpts.boundaryStepMinRel*gammaSuccess);
        boundaryStep = max(boundaryStepMin, ...
            searchOpts.boundaryStepRel*gammaSuccess);

        for ib = 2:searchOpts.maxBoundaryTrials
            gammaBoundary = gammaSuccess + boundaryStep;
            warmBoundary = make_warm_init(initState,solSuccess);

            fprintf(['    [%s-range] full-budget boundary %02d: ', ...
                'Gamma=%.6e (%.3f dB), step=%.6e\n'], ...
                rangeTag,ib-1,gammaBoundary,10*log10(gammaBoundary),boundaryStep);

            [solBoundary,histBoundary] = run_trial(params,ch,alg, ...
                warmBoundary,gammaBoundary,fullBudget);
            boundaryTrialsUsed = boundaryTrialsUsed + 1;
            [trace,trialCounter] = add_record(trace,trialCounter, ...
                'boundary-extend',gammaBoundary,fullBudget, ...
                solBoundary,boundaryStep);

            if solBoundary.success
                gammaSuccess = gammaBoundary;
                solSuccess = solBoundary;
                histSuccess = histBoundary;

                % A previously failed short-budget point below the newly
                % confirmed success is no longer a valid upper boundary.
                if isfinite(gammaFail) && gammaSuccess >= gammaFail*(1-1e-10)
                    gammaFail = NaN;
                    upperBoundFound = false;
                end
                continue;
            end

            boundaryFailedUpper = gammaBoundary;
            gammaFail = gammaBoundary;
            upperBoundFound = true;

            nextStep = boundaryStep*searchOpts.boundaryStepShrink;
            if nextStep < boundaryStepMin*(1-1e-12)
                boundaryFinalStep = boundaryStep;
                break;
            end

            boundaryStep = nextStep;
            boundaryFinalStep = boundaryStep;
            fprintf(['    [%s-range] full-budget boundary failed; ', ...
                'reducing step to %.6e.\n'],rangeTag,boundaryStep);
        end

        if isnan(boundaryFinalStep)
            boundaryFinalStep = boundaryStep;
        end
    else
        fprintf(['    [%s-range] full-budget confirmation did not ', ...
            'produce a retained feasible pair; keeping the earlier ', ...
            'verified boundary candidate.\n'],rangeTag);
    end
end

gammaMaxFound = gammaSuccess;
solAtGammaMax = solSuccess;

if upperBoundFound && isfinite(gammaFail)
    finalRelativeGap = (gammaFail-gammaSuccess)/max(gammaSuccess,realmin);
else
    finalRelativeGap = NaN;
end

searchInfo.success = true;
searchInfo.status = 'A feasible SIM SINR range was found.';
searchInfo.searchMode = searchMode;
searchInfo.gammaRequestedStart = gammaStart;
searchInfo.gammaStartFeasible = gammaStartFeasible;
searchInfo.gammaMaxFound = gammaMaxFound;
searchInfo.gammaFailedUpper = gammaFail;
searchInfo.upperBoundFound = upperBoundFound;
searchInfo.finalRelativeGap = finalRelativeGap;
searchInfo.recoveryTrialsUsed = recoveryUsed;
searchInfo.expandTrialsUsed = expandUsed;
searchInfo.confirmTrialsUsed = confirmUsed;
searchInfo.bridgeTrialsUsed = bridgeUsed;
searchInfo.stepReductionsUsed = stepReductionUsed;
searchInfo.refineTrialsUsed = refineUsed;
searchInfo.gammaBeforeBoundaryConfirm = gammaBeforeBoundaryConfirm;
searchInfo.boundaryTrialsUsed = boundaryTrialsUsed;
searchInfo.boundaryFailedUpper = boundaryFailedUpper;
searchInfo.boundaryFinalStep = boundaryFinalStep;
searchInfo.fullBudgetConfirmed = fullBudgetConfirmed;
searchInfo.options = searchOpts;
searchInfo.trace = trace;
searchInfo.histAtGammaMax = histSuccess;
searchInfo.bestTheta = solSuccess.theta;
if isfield(solSuccess,'W')
    searchInfo.bestW = solSuccess.W;
else
    searchInfo.bestW = [];
end

fprintf(['    [%s-range] largest verified Gamma=%.6e (%.3f dB), ', ...
    'start recovered to %.6e, upperBoundFound=%d'],rangeTag,gammaMaxFound, ...
    10*log10(gammaMaxFound),gammaStartFeasible,upperBoundFound);
if upperBoundFound
    fprintf(', final bracket gap=%.3f%%',100*finalRelativeGap);
end
fprintf(', fullBudgetConfirmed=%d, boundaryTrials=%d\n', ...
    fullBudgetConfirmed,boundaryTrialsUsed);
end

function [sol,hist] = run_trial(params,ch,alg,initState,gammaValue,maxOuter)
algTrial = alg;
algTrial.maxOuter = maxOuter;
gammaTarget = gammaValue*ones(params.K,1);
[sol,hist] = SIM_run_CRB_SINR_AO_solver( ...
    params,ch,algTrial,initState,gammaTarget);
validate_solver_contract(sol);
end

function validate_solver_contract(sol)
requiredFields = {'success','status','bestIteration','bestStage', ...
    'failureStage','firstSDRFeasible','firstSDRStatus'};
for ii = 1:numel(requiredFields)
    name = requiredFields{ii};
    if ~isfield(sol,name)
        error('SIM_search_max_common_sinr_AO:OldAOSolverDetected', ...
            ['The active SIM_run_CRB_SINR_AO_solver.m does not return the ', ...
             'required field "%s". Replace it with the revised AO ', ...
             'solver used by the current CRB-SINR scripts.'],name);
    end
end
end

function warm = make_warm_init(baseInit,sol)
warm = baseInit;
if isfield(sol,'theta') && ~isempty(sol.theta)
    warm.theta0 = sol.theta;
end
if isfield(sol,'W') && ~isempty(sol.W)
    warm.W0 = sol.W;
end
if isfield(sol,'QcommSDR') && ~isempty(sol.QcommSDR), warm.Qcomm0=sol.QcommSDR; end
if isfield(sol,'RsSDR') && ~isempty(sol.RsSDR), warm.Rs0=sol.RsSDR; end
if isfield(sol,'S') && ~isempty(sol.S)
    warm.S = sol.S;
end
if isfield(sol,'sigma_c2') && ~isempty(sol.sigma_c2)
    warm.sigma_c2 = sol.sigma_c2;
end
if isfield(sol,'sigma_s2') && ~isempty(sol.sigma_s2)
    warm.sigma_s2 = sol.sigma_s2;
end
end

function tf = is_first_sdr_failure(sol)
tf = isstruct(sol) && isfield(sol,'firstSDRFeasible') ...
    && isscalar(sol.firstSDRFeasible) && ~logical(sol.firstSDRFeasible);
end

function [trace,count] = add_record(trace,count,stage,gammaValue,budget,sol,step)
count = count + 1;
trace(count).stage = stage;
trace(count).gamma = gammaValue;
trace(count).gammaDB = 10*log10(max(gammaValue,realmin));
trace(count).step = step;
trace(count).maxOuter = budget;
trace(count).success = isfield(sol,'success') && isscalar(sol.success) ...
    && logical(sol.success);
trace(count).status = status_text(sol,'status');
trace(count).failureStage = status_text(sol,'failureStage');
trace(count).firstSDRFeasible = false;
if isfield(sol,'firstSDRFeasible') && isscalar(sol.firstSDRFeasible)
    trace(count).firstSDRFeasible = logical(sol.firstSDRFeasible);
end
trace(count).firstSDRStatus = status_text(sol,'firstSDRStatus');
trace(count).bestIteration = NaN;
trace(count).bestStage = 'none';
trace(count).CRB = NaN;
trace(count).minSINR = NaN;
if isfield(sol,'bestIteration') && isscalar(sol.bestIteration)
    trace(count).bestIteration = sol.bestIteration;
end
if isfield(sol,'bestStage') && ~isempty(sol.bestStage)
    trace(count).bestStage = char(string(sol.bestStage));
end
if trace(count).success && isfield(sol,'metrics')
    if isfield(sol.metrics,'CRB') && isscalar(sol.metrics.CRB)
        trace(count).CRB = sol.metrics.CRB;
    end
    if isfield(sol.metrics,'sinr') && ~isempty(sol.metrics.sinr)
        trace(count).minSINR = min(sol.metrics.sinr);
    end
end
end

function textValue = status_text(s,fieldName)
textValue = 'not-reported';
if isstruct(s) && isfield(s,fieldName) && ~isempty(s.(fieldName))
    textValue = char(string(s.(fieldName)));
end
end

function searchInfo = make_failure_search_info(opts,trace,gammaStart,lastHist,status)
searchInfo.success = false;
searchInfo.status = status;
searchInfo.searchMode = lower(char(string(opts.searchMode)));
searchInfo.gammaRequestedStart = gammaStart;
searchInfo.gammaStartFeasible = NaN;
searchInfo.gammaMaxFound = NaN;
searchInfo.gammaFailedUpper = NaN;
searchInfo.upperBoundFound = false;
searchInfo.finalRelativeGap = NaN;
searchInfo.recoveryTrialsUsed = NaN;
searchInfo.expandTrialsUsed = 0;
searchInfo.confirmTrialsUsed = 0;
searchInfo.bridgeTrialsUsed = 0;
searchInfo.stepReductionsUsed = 0;
searchInfo.refineTrialsUsed = 0;
searchInfo.gammaBeforeBoundaryConfirm = NaN;
searchInfo.boundaryTrialsUsed = 0;
searchInfo.boundaryFailedUpper = NaN;
searchInfo.boundaryFinalStep = NaN;
searchInfo.fullBudgetConfirmed = false;
searchInfo.options = opts;
searchInfo.trace = trace;
searchInfo.histAtGammaMax = lastHist;
searchInfo.bestTheta = [];
searchInfo.bestW = [];
end

function sol = make_empty_failure_solution(gammaValue)
sol.success = false;
sol.status = 'No trial has been run.';
sol.bestIteration = NaN;
sol.bestStage = 'none';
sol.failureStage = 'not-run';
sol.firstSDRFeasible = false;
sol.firstSDRStatus = 'not-run';
sol.gammaTarget = gammaValue;
sol.theta = [];
sol.W = [];
sol.metrics = struct();
end

function opts = fill_defaults(opts)
defaults.searchMode = 'multiplicative';
defaults.growthFactor = 1.2;
defaults.additiveStep = 0.1;
defaults.additiveStepMin = 1e-3;
defaults.stepShrink = 0.5;
defaults.maxStepHalving = 10;
defaults.maxOuterGrowth = 30;
defaults.maxOuterConfirm = 80;
defaults.maxOuterRefine = 50;
defaults.maxExpand = 100;
defaults.maxBridge = 3;
defaults.maxRefine = 8;
defaults.rangeRelTol = 0.01;
defaults.enableRecovery = false;
defaults.recoveryFactor = 0.9;
defaults.maxRecovery = 30;
defaults.recoveryMinGamma = 1e-8;
defaults.enableFullBudgetBoundary = true;
defaults.maxBoundaryTrials = 12;
defaults.boundaryStepRel = 0.01;
defaults.boundaryStepMinRel = 1e-4;
defaults.boundaryStepMinAbs = 1e-4;
defaults.boundaryStepShrink = 0.5;

names = fieldnames(defaults);
for ii = 1:numel(names)
    name = names{ii};
    if ~isfield(opts,name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
end
