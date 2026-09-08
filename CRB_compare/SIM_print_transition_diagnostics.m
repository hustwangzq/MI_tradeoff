function SIM_print_transition_diagnostics(mode,diagnosticsOverride, ...
    resultOverride,reportOverride)
%SIM_PRINT_TRANSITION_DIAGNOSTICS Print the saved dual-state/refinement trace.
%
% Command-line use:
%   SIM_print_transition_diagnostics
%   SIM_print_transition_diagnostics('all')
%   SIM_print_transition_diagnostics('abnormal')
%   SIM_print_transition_diagnostics('summary')
%
% The default is 'all'.  The same text printed in the command window is
% also written to SIM_CRB_SINR_diagnostics_report.txt in the current MATLAB working directory.

if nargin < 1 || isempty(mode)
    mode = 'all';
end
mode = lower(char(string(mode)));
if ~ismember(mode,{'all','abnormal','summary'})
    error('SIM_print_transition_diagnostics:BadMode', ...
        'Mode must be ''all'', ''abnormal'', or ''summary''.');
end

thisFile = mfilename('fullpath');
if isempty(thisFile)
    packageDir = pwd;
else
    packageDir = fileparts(thisFile);
end
resultsDir = fullfile(packageDir,'results');
diagnosticsFile = fullfile(resultsDir,'SIM_CRB_SINR_diagnostics.mat');
resultFile = fullfile(resultsDir,'CRB_SINR_SIM_results.mat');
reportFile = fullfile(pwd,'SIM_CRB_SINR_diagnostics_report.txt');
if nargin>=2&&~isempty(diagnosticsOverride)
    diagnosticsFile=char(string(diagnosticsOverride));
end
if nargin>=3&&~isempty(resultOverride)
    resultFile=char(string(resultOverride));
end
if nargin>=4&&~isempty(reportOverride)
    reportFile=char(string(reportOverride));
end

if exist(diagnosticsFile,'file') ~= 2
    error('SIM_print_transition_diagnostics:MissingDiagnostics', ...
        ['Diagnostics file not found:\n%s\nRun SIM_CRB_SINR first.'], ...
        diagnosticsFile);
end

loaded = load(diagnosticsFile,'transitionDiagnostics','scriptVersion', ...
    'dualState','refine','diagnosticOptions');
if ~isfield(loaded,'transitionDiagnostics') || ...
        ~istable(loaded.transitionDiagnostics)
    error('SIM_print_transition_diagnostics:InvalidDiagnostics', ...
        'The diagnostics MAT file does not contain a valid table.');
end
T = loaded.transitionDiagnostics;

fid = fopen(reportFile,'w');
if fid < 0
    error('SIM_print_transition_diagnostics:OpenReportFailed', ...
        'Cannot open report file for writing: %s',reportFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

emit(fid,'============================================================\n');
emit(fid,'SIM protected-continuation / robust-CVX diagnostics\n');
emit(fid,'============================================================\n');
if isfield(loaded,'scriptVersion')
    emit(fid,'Script version : %s\n',char(loaded.scriptVersion));
end
emit(fid,'Diagnostics extension: FIM/continuity/solver tracing enabled\n');
emit(fid,'Required continuity fields: deltaTheta deltaW deltaP branchSwitchFlag\n');
emit(fid,'Rows saved     : %d\n',height(T));
emit(fid,'Diagnostics MAT: %s\n',diagnosticsFile);
emit(fid,'Text report    : %s\n\n',reportFile);
if isfield(loaded,'diagnosticOptions')&& ...
        isfield(loaded.diagnosticOptions,'traceTextFile')
    emit(fid,'Per-point root-cause TXT: %s\n\n', ...
        char(loaded.diagnosticOptions.traceTextFile));
end

if isempty(T)
    emit(fid,'No diagnostic rows were saved.\n');
    return;
end

layers = unique(T.Layer(isfinite(T.Layer))).';
for L = layers
    TL = T(T.Layer == L,:);
    bridgeRows = TL.BridgeAttempted;
    handoffText = cell_text_column(TL.HandoffSource);
    bridgeAccepted = sum(startsWith(handoffText,'bridge-'));
    reportFallback = sum(strcmp(handoffText,'report-fallback'));
    robustRows = 0;
    fallbackRows = 0;
    if ismember('CVXRobustUsed',TL.Properties.VariableNames)
        robustRows = sum(TL.CVXRobustUsed);
    end
    if ismember('CVXAttempt',TL.Properties.VariableNames)
        fallbackRows = sum(isfinite(TL.CVXAttempt) & TL.CVXAttempt > 1);
    end
    refinementCount = sum(TL.RefinementInserted);
    severityText = cell_text_column(TL.Severity);
    hardCount = sum(strcmp(severityText,'hard'));
    softCount = sum(strcmp(severityText,'soft'));
    failureCount = sum(strcmp(severityText,'failure'));
    leftBendCount = sum(isfinite(TL.DeltaCRBDB) & TL.DeltaCRBDB < -0.10);
    positiveJumpCount = sum(isfinite(TL.DeltaCRBDB) & TL.DeltaCRBDB > 1.0);
    maxPositiveJump = safe_max(TL.DeltaCRBDB(TL.DeltaCRBDB > 0));
    maxLeftBend = safe_max(-TL.DeltaCRBDB(TL.DeltaCRBDB < 0));

    emit(fid,'======================== SIM L=%d ========================\n',L);
    emit(fid,'Target-solve attempts       : %d\n',height(TL));
    emit(fid,'Rows with bridge attempted  : %d\n',sum(bridgeRows));
    emit(fid,'Bridge handoff accepted     : %d\n',bridgeAccepted);
    emit(fid,'Report fallback after protected handoff: %d\n',reportFallback);
    emit(fid,'Rows using robust CVX       : %d\n',robustRows);
    emit(fid,'Rows needing CVX fallback   : %d\n',fallbackRows);
    emit(fid,'Adaptive refinements inserted: %d\n',refinementCount);
    emit(fid,'Soft / hard / failure rows  : %d / %d / %d\n', ...
        softCount,hardCount,failureCount);
    emit(fid,'Positive jumps > 1 dB       : %d\n',positiveJumpCount);
    emit(fid,'Left bends > 0.1 dB         : %d\n',leftBendCount);
    emit(fid,'Largest positive jump       : %.3f dB\n',maxPositiveJump);
    emit(fid,'Largest apparent left bend  : %.3f dB\n',maxLeftBend);
    print_dual_path_layer_summary(fid,TL);
    emit(fid,'\n');
end

maxAbnormalDetails=12;
if isfield(loaded,'diagnosticOptions')&& ...
        isfield(loaded.diagnosticOptions,'maxAbnormalDetailsPerLayer')
    maxAbnormalDetails=loaded.diagnosticOptions.maxAbnormalDetailsPerLayer;
end
print_bounded_root_cause_details(fid,T,maxAbnormalDetails);

% Keep the new outer-margin diagnostics compact: one line per layer only.
print_outer_margin_summary(fid,resultFile);

% Keep deferred rank-one diagnostics deliberately compact: one summary line
% per layer, followed only by a bounded number of failed-point details.
print_deferred_recovery_summary(fid,resultFile,5);

if strcmp(mode,'summary')
    emit(fid,'Summary-only mode selected.\n');
    return;
end

emit(fid,'====================== Detailed trace ======================\n');
for ii = 1:height(T)
    isAbnormal = T.RefinementInserted(ii) || T.BranchSwitchTried(ii) || ...
        ~strcmp(cell_text(T.Severity,ii),'none') || ...
        strcmp(cell_text(T.HandoffSource,ii),'report-fallback') || ...
        (isfinite(T.DeltaCRBDB(ii)) && abs(T.DeltaCRBDB(ii)) > 1.0);
    if strcmp(mode,'abnormal') && ~isAbnormal
        continue;
    end

    emit(fid,'\n[%04d] SIM L=%d | coarse=%d | %s | depth=%d\n', ...
        ii,T.Layer(ii),T.CoarseIndex(ii), ...
        cell_text(T.TargetKind,ii),T.RefinementDepth(ii));
    emit(fid,'  target path : %.3f -> %.3f dB  (step %+0.3f dB)\n', ...
        T.PrevTargetDB(ii),T.TargetDB(ii),T.StepDB(ii));
    emit(fid,'  grid flags  : forced=%d, inserted=%d, dualWindow=%d\n', ...
        T.IsForced(ii),T.IsInserted(ii),T.DualStateWindow(ii));

    emit(fid,'  handoff     : %s\n',cell_text(T.HandoffSource,ii));
    if T.BridgeAttempted(ii)
        emit(fid,'  protected   : attempts=%d, progress=%.2f, status=%s\n', ...
            T.BridgeNumAttempts(ii),T.BridgeEtaUsed(ii), ...
            cell_text(T.BridgeStatus,ii));
        if startsWith(cell_text(T.HandoffSource,ii),'bridge-')
            emit(fid,'                accepted target=%.3f dB, CRB=%.3f dB, ', ...
                T.BridgeTargetDB(ii),T.BridgeCRBDB(ii));
            emit(fid,'SINR=%.3f dB, relEig(A/B)=%.2e/%.2e\n', ...
                T.BridgeMinSINRDB(ii),T.BridgeRelEigA(ii),T.BridgeRelEigB(ii));
        else
            emit(fid,'                no protected step accepted; report state was handed off.\n');
        end
        emit(fid,'  protected tries: %s\n',cell_text(T.BridgeAttemptSummary,ii));
    else
        emit(fid,'  protected   : not attempted\n');
    end

    if T.ReportSuccess(ii)
        emit(fid,'  report try  : CRB=%.3f dB, minSINR=%.3f dB, ', ...
            T.ReportCRBDB(ii),T.ReportMinSINRDB(ii));
        emit(fid,'relEig(A/B)=%.2e/%.2e\n', ...
            T.ReportRelEigA(ii),T.ReportRelEigB(ii));
    else
        emit(fid,'  report try  : failed / no retained complete candidate\n');
    end
    if ismember('CVXStatus',T.Properties.VariableNames)
        emit(fid,['  CVX         : robust=%d, attempt=%.0f, status=%s, ', ...
            'tol=%.2e, precision=%s, scale=%.2f, reg=%.1e\n'], ...
            T.CVXRobustUsed(ii),T.CVXAttempt(ii),cell_text(T.CVXStatus,ii), ...
            T.CVXSlvtol(ii),cell_text(T.CVXPrecision,ii), ...
            T.CVXSchurScale(ii),T.CVXRegRel(ii));
    end
    emit(fid,'  transition  : dCRB=%+.3f dB, severity=%s\n', ...
        T.DeltaCRBDB(ii),cell_text(T.Severity,ii));

    % Extra diagnostics (printed only when saved in the table).
    if ismember('lambdaMinA',T.Properties.VariableNames)
        emit(fid,'  FIM         : lambdaA=%.3e lambdaB=%.3e condA=%.3e condB=%.3e rankA=%g rankB=%g\n', ...
            T.lambdaMinA(ii),T.lambdaMinB(ii),T.condA(ii),T.condB(ii),T.rankA(ii),T.rankB(ii));
    end
    if ismember('deltaTheta',T.Properties.VariableNames)
        emit(fid,'  continuity  : deltaTheta=%.3e deltaW=%.3e deltaP=%.3e\n', ...
            T.deltaTheta(ii),T.deltaW(ii),T.deltaP(ii));
    end

    if T.RefinementInserted(ii)
        emit(fid,'  action      : INSERT midpoint %.3f dB and retry target\n', ...
            T.InsertedMidDB(ii));
    else
        emit(fid,'  action      : no new midpoint inserted\n');
    end
    if T.BranchSwitchTried(ii)
        emit(fid,'  branch      : branch recovery was attempted\n');
    end
    emit(fid,'  source      : %s\n',cell_text(T.SolutionSource,ii));
    emit(fid,'  recovery    : %s, covarianceGap=%.3e\n', ...
        cell_text(T.RecoveryStatus,ii),T.RecoveryCovarianceGapRel(ii));
    reasonText = cell_text(T.RecoveryReason,ii);
    if ~isempty(reasonText)
        emit(fid,'  reason      : %s\n',reasonText);
    end
end

emit(fid,'\n============================================================\n');
emit(fid,'Finished. The same report was saved to:\n%s\n',reportFile);
end

function print_dual_path_layer_summary(fid,T)
required={'ContinuationPath','CarrierPath','CarrierQuality', ...
    'RootCause','CarrierActualNextMargin','CarrierTargetNextMargin', ...
    'TotalSolveSeconds'};
if ~all(ismember(required,T.Properties.VariableNames)),return;end
primaryAttempt=0;primarySuccess=0;primaryFirstFail=0;
carrierAttempt=0;carrierSuccess=0;carrierFirstFail=0;
innerBlockedSINR=0;innerBlockedObjective=0;innerBlockedCondition=0;
for ii=1:height(T)
    primary=T.ContinuationPath{ii};
    carrier=T.CarrierPath{ii};
    primaryAttempt=primaryAttempt+logical(primary.attempted);
    primarySuccess=primarySuccess+logical(primary.success);
    primaryFirstFail=primaryFirstFail+ ...
        (logical(primary.attempted)&&~logical(primary.firstSDRFeasible));
    carrierAttempt=carrierAttempt+logical(carrier.attempted);
    carrierSuccess=carrierSuccess+logical(carrier.success);
    carrierFirstFail=carrierFirstFail+ ...
        (logical(carrier.attempted)&&~logical(carrier.firstSDRFeasible));
    innerBlockedSINR=innerBlockedSINR+ ...
        (primary.phaseAccepted==0&&primary.phaseFeasibilityRejects>0);
    innerBlockedObjective=innerBlockedObjective+ ...
        (primary.phaseAccepted==0&&primary.phaseObjectiveRejects>0);
    innerBlockedCondition=innerBlockedCondition+ ...
        (primary.phaseAccepted==0&&primary.phaseConditionRejects>0);
end
quality=cell_text_column(T.CarrierQuality);
carrierSelected=sum(strcmp(cell_text_column(T.SolutionSource), ...
    'persistent-commCarrier'));
reserveAvailable=isfinite(T.CarrierActualNextMargin)& ...
    isfinite(T.CarrierTargetNextMargin);
reserveShortfall=sum(reserveAvailable& ...
    T.CarrierActualNextMargin<T.CarrierTargetNextMargin);
elapsed=T.TotalSolveSeconds(isfinite(T.TotalSolveSeconds));
emit(fid,'bestCRB AO attempted / solved / first-SDR-failed: %d / %d / %d\n', ...
    primaryAttempt,primarySuccess,primaryFirstFail);
emit(fid,'Carrier AO attempted / solved / first-SDR-failed: %d / %d / %d\n', ...
    carrierAttempt,carrierSuccess,carrierFirstFail);
emit(fid,'Carrier won / inherited / reserve-shortfall: %d / %d / %d\n', ...
    carrierSelected,sum(T.CarrierUsedPrevious),reserveShortfall);
emit(fid,'Carrier quality protected / relaxed / fallback: %d / %d / %d\n', ...
    sum(quality=="protected"),sum(quality=="relaxed"), ...
    sum(quality=="fallback"));
emit(fid,'Inner-phase blocked SINR / CRB / conditioning: %d / %d / %d\n', ...
    innerBlockedSINR,innerBlockedObjective,innerBlockedCondition);
if ~isempty(elapsed)
    emit(fid,'Measured AO+carrier time total / median / max: %.1f min / %.1f s / %.1f s\n', ...
        sum(elapsed)/60,median(elapsed),max(elapsed));
end
causes=cell_text_column(T.RootCause);
interesting=causes~="normal-continuation"& ...
    causes~="carrier-improved-bestCRB";
if any(interesting)
    selectedCauses=causes(interesting);
    uniqueCauses=unique(selectedCauses,'stable');
    for jj=1:numel(uniqueCauses)
        emit(fid,'  root-cause %-42s %d\n', ...
            char(uniqueCauses(jj)),sum(selectedCauses==uniqueCauses(jj)));
    end
end
end

function print_bounded_root_cause_details(fid,T,maxPerLayer)
if ~all(ismember({'RootCause','ContinuationPath','CarrierPath', ...
        'CarrierActualNextMargin','CarrierTargetNextMargin', ...
        'CarrierQuality','CarrierCRBDeltaDB'},T.Properties.VariableNames))
    return;
end
emit(fid,'================ Bounded jump/root-cause details ===============\n');
layers=unique(T.Layer(isfinite(T.Layer))).';
for L=layers
    layerRows=find(T.Layer==L);
    severities=cell_text_column(T.Severity(layerRows));
    causes=cell_text_column(T.RootCause(layerRows));
    scores=zeros(numel(layerRows),1);
    delta=T.DeltaCRBDB(layerRows);
    finiteDelta=isfinite(delta);
    scores(finiteDelta)=abs(delta(finiteDelta));
    scores= scores+100*(severities=="failure")+ ...
        40*(severities=="hard")+20*(severities=="soft")+ ...
        10*(causes=="carrier-feasible-only-with-large-CRB-loss");
    interesting=severities~="none"| ...
        (finiteDelta&abs(delta)>0.10)| ...
        strcmp(cell_text_column(T.CarrierQuality(layerRows)),"fallback");
    ranked=find(interesting);
    if isempty(ranked)
        emit(fid,'L=%d: no abnormal dual-path transition recorded.\n',L);
        continue;
    end
    [~,order]=sort(scores(ranked),'descend');
    ranked=ranked(order);
    shown=min(numel(ranked),maxPerLayer);
    emit(fid,'L=%d: showing %d of %d abnormal attempts, ranked by severity.\n', ...
        L,shown,numel(ranked));
    for jj=1:shown
        ii=layerRows(ranked(jj));
        primary=T.ContinuationPath{ii};
        carrier=T.CarrierPath{ii};
        emit(fid,['  target %.5f dB | dCRB=%+.3f dB | %s | root=%s | ', ...
            'reserve=%+.3e/%.3e | carrierDelta=%+.2f dB [%s]\n'], ...
            T.TargetDB(ii),T.DeltaCRBDB(ii),T.Severity{ii}, ...
            T.RootCause{ii},T.CarrierActualNextMargin(ii), ...
            T.CarrierTargetNextMargin(ii),T.CarrierCRBDeltaDB(ii), ...
            T.CarrierQuality{ii});
        emit(fid,['    bestCRB success=%d firstSDR=%d CVX=%s CRB=%.3f ', ...
            'outer=%d PhiAcc=%d inner=%d reject(SINR/CRB/FIM)=%d/%d/%d\n'], ...
            primary.success,primary.firstSDRFeasible,primary.firstCVXStatus, ...
            primary.physicalCRBDB,primary.outerIterations, ...
            primary.phaseAccepted,primary.innerPhaseSteps, ...
            primary.phaseFeasibilityRejects,primary.phaseObjectiveRejects, ...
            primary.phaseConditionRejects);
        emit(fid,['    carrier success=%d firstSDR=%d CVX=%s CRB=%.3f ', ...
            'outer=%d PhiAcc=%d inner=%d reject(SINR/CRB/FIM)=%d/%d/%d\n'], ...
            carrier.success,carrier.firstSDRFeasible,carrier.firstCVXStatus, ...
            carrier.physicalCRBDB,carrier.outerIterations, ...
            carrier.phaseAccepted,carrier.innerPhaseSteps, ...
            carrier.phaseFeasibilityRejects,carrier.phaseObjectiveRejects, ...
            carrier.phaseConditionRejects);
    end
    if numel(ranked)>shown
        emit(fid,'  %d additional attempts remain in the MAT/per-point TXT.\n', ...
            numel(ranked)-shown);
    end
end
emit(fid,'\n');
end

function emit(fid,varargin)
fprintf(varargin{:});
fprintf(fid,varargin{:});
end

function out = cell_text_column(v)
out = strings(numel(v),1);
for ii = 1:numel(v)
    out(ii) = string(cell_text(v,ii));
end
end

function text = cell_text(v,idx)
value = v{idx};
if iscell(value) && isscalar(value)
    value = value{1};
end
if isstring(value)
    text = char(value);
elseif ischar(value)
    text = value;
else
    text = char(string(value));
end
end

function value = safe_max(x)
if isempty(x)
    value = 0;
else
    value = max(x);
end
end

function print_deferred_recovery_summary(fid,resultFile,maxFailureDetails)
emit(fid,'================ Deferred rank-one recovery ================\n');
if exist(resultFile,'file') ~= 2
    emit(fid,'Result MAT unavailable; deferred recovery summary skipped.\n\n');
    return;
end

data = load(resultFile,'simCRBResults');
if ~isfield(data,'simCRBResults') || isempty(data.simCRBResults)
    emit(fid,'No saved SIM results; deferred recovery summary skipped.\n\n');
    return;
end

failures = struct('L',{},'gammaDB',{},'status',{});
for idx = 1:numel(data.simCRBResults)
    entry = data.simCRBResults(idx);
    tr = entry.trade;
    total = 0; okCount = 0; failCount = 0;
    covGap = []; powerGap = []; desiredGap = []; leakGap = [];
    for ip = 1:numel(tr)
        if ~isfield(tr(ip),'success') || ~logical(tr(ip).success) || ...
                ~isfield(tr(ip),'rankOneRecovery') || ...
                isempty(tr(ip).rankOneRecovery)
            continue;
        end
        total = total + 1;
        rec = tr(ip).rankOneRecovery;
        valid = isfield(rec,'valid') && isscalar(rec.valid) && logical(rec.valid);
        if valid
            okCount = okCount + 1;
        else
            failCount = failCount + 1;
            item.L = entry.L;
            item.gammaDB = scalar_or_nan(tr(ip),'targetGammaDB');
            item.status = text_or_default(rec,'status','unknown');
            failures(end+1) = item; %#ok<AGROW>
        end
        covGap(end+1) = scalar_or_nan(rec,'covarianceGapRel'); %#ok<AGROW>
        powerGap(end+1) = scalar_or_nan(rec,'powerGapRel'); %#ok<AGROW>
        desiredGap(end+1) = scalar_or_nan(rec,'maxDesiredSignalGapRel'); %#ok<AGROW>
        leakGap(end+1) = scalar_or_nan(rec,'maxSelfResidualLeakRel'); %#ok<AGROW>
    end
    emit(fid,['L=%d: recovered=%d, valid=%d, failed=%d, ', ...
        'maxGap(cov/power/desired/leak)=%.2e/%.2e/%.2e/%.2e\n'], ...
        entry.L,total,okCount,failCount,finite_max(covGap), ...
        finite_max(powerGap),finite_max(desiredGap),finite_max(leakGap));
end

shown = min(numel(failures),maxFailureDetails);
for ii = 1:shown
    emit(fid,'  FAILED L=%d, target=%.3f dB: %s\n', ...
        failures(ii).L,failures(ii).gammaDB,failures(ii).status);
end
if numel(failures) > shown
    emit(fid,'  ... %d additional failed recovery point(s) omitted.\n', ...
        numel(failures)-shown);
end
emit(fid,'\n');
end

function print_outer_margin_summary(fid,resultFile)
emit(fid,'================ Outer-phase SINR reserve ===================\n');
if exist(resultFile,'file') ~= 2
    emit(fid,'Result MAT unavailable; outer-margin summary skipped.\n\n');
    return;
end
data = load(resultFile,'simCRBResults');
if ~isfield(data,'simCRBResults') || isempty(data.simCRBResults)
    emit(fid,'No saved SIM results; outer-margin summary skipped.\n\n');
    return;
end
for idx = 1:numel(data.simCRBResults)
    entry = data.simCRBResults(idx);
    preApplied = 0; outerSteps = 0; innerSteps = 0; points = 0;
    for ip = 1:numel(entry.trade)
        tr = entry.trade(ip);
        if ~isfield(tr,'success') || ~logical(tr.success) || ...
                ~isfield(tr,'hist') || ~isstruct(tr.hist)
            continue;
        end
        points = points+1;
        h = tr.hist;
        if isfield(h,'preOuterMarginApplied') && ...
                isscalar(h.preOuterMarginApplied)
            preApplied = preApplied+logical(h.preOuterMarginApplied);
        end
        if isfield(h,'outerMarginAcceptedSteps')
            outerSteps = outerSteps+sum(h.outerMarginAcceptedSteps(:),'omitnan');
        end
        if isfield(h,'phaseInner')
            innerSteps = innerSteps+sum(h.phaseInner(:),'omitnan');
        end
    end
    emit(fid,['L=%d: successfulPoints=%d, preTargetOuterApplied=%d, ', ...
        'outerMarginSteps=%d, innerCRBSteps=%d\n'], ...
        entry.L,points,preApplied,round(outerSteps),round(innerSteps));
end
emit(fid,'\n');
end

function value = scalar_or_nan(s,name)
value = NaN;
if isstruct(s) && isfield(s,name) && isscalar(s.(name)) && ...
        isnumeric(s.(name)) && isfinite(s.(name))
    value = s.(name);
end
end

function value = finite_max(x)
x = x(isfinite(x));
if isempty(x), value = NaN; else, value = max(x); end
end

function value = text_or_default(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = char(string(s.(name)));
end
end
