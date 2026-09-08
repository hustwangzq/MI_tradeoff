function [gammaMax,sol,info] = SIM_max_common_sinr_fixed_phase_SDR(params,ch,P,opts)
%SIM_MAX_COMMON_SINR_FIXED_PHASE_SDR Pure-communication common-SINR bound.
%
% For a FIXED SIM phase matrix P, this function maximizes the common SINR
% without any sensing/CRB constraint. It is used only to estimate the
% communication-side ceiling. Across multiple phase initializations, the
% largest fixed-P result is a reproducible lower bound on the joint
% pure-communication optimum.

if nargin < 4 || isempty(opts)
    opts = struct('gammaSeed',1);
end

K = params.K;
Nt = params.Nt;
P0 = params.P0;
maxExpand = max(5,round(get_option(opts,'maxExpand',35)));
maxBisection = max(10,round(get_option(opts,'maxBisection',32)));
gammaSeed = get_option(opts,'gammaSeed',1);
gammaSeed = max(gammaSeed,1e-8);
upperCapDB = get_option(opts,'upperCapDB',45);
upperCap = 10^(upperCapDB/10);
relTol = get_option(opts,'relativeTolerance',2e-4);
cvxQuiet = logical(get_option(opts,'cvxQuiet',true));
precisionMode = lower(char(string(get_option(opts,'cvxPrecisionMode','default'))));
if ~ismember(precisionMode,{'default','high','best'})
    precisionMode = 'default';
end

% Fixed effective user channels.
gUsers = zeros(Nt,K);
Hscaled = zeros(Nt,Nt,K);
noiseScaled = zeros(K,1);
for k = 1:K
    gk = ch.G' * P' * ch.hUsers{k};
    gk = gk(:);
    gUsers(:,k) = gk;
    Hk = (gk*gk' + (gk*gk')')/2;
    userScale = max(P0*real(trace(Hk))+ch.sigma_c2,realmin);
    Hscaled(:,:,k) = (P0/userScale)*Hk;
    noiseScaled(k) = ch.sigma_c2/userScale;
end

low = 0;
high = gammaSeed;
bestQ = zeros(Nt,Nt,K);
lastStatus = 'not-run';
numFeasibilityChecks = 0;

% First locate one feasible positive SINR.
[seedOK,Qseed,statusSeed] = fixed_phase_feasibility( ...
    gammaSeed,Nt,K,P0,Hscaled,noiseScaled,gUsers,ch.sigma_c2,opts, ...
    cvxQuiet,precisionMode);
numFeasibilityChecks = numFeasibilityChecks + 1;
lastStatus = statusSeed;
if seedOK
    low = gammaSeed;
    bestQ = Qseed;
else
    trial = gammaSeed;
    for ishrink = 1:30
        trial = 0.5*trial;
        [ok,Qtmp,statusTmp] = fixed_phase_feasibility( ...
            trial,Nt,K,P0,Hscaled,noiseScaled,gUsers,ch.sigma_c2,opts, ...
            cvxQuiet,precisionMode);
        numFeasibilityChecks = numFeasibilityChecks + 1;
        lastStatus = statusTmp;
        if ok
            low = trial;
            bestQ = Qtmp;
            break;
        end
    end
end

if low == 0
    gammaMax = NaN;
    sol = make_solution(false,NaN,bestQ,zeros(K,1),'no-feasible-positive-SINR');
    info = make_info(false,numFeasibilityChecks,lastStatus,NaN,NaN,precisionMode);
    return;
end

% Expand upward from the known feasible lower bound until an infeasible
% upper bracket is found or the configured physical cap is reached.
high = min(max(gammaSeed,2*low),upperCap);
highFeasible = false;
for ie = 1:maxExpand
    [ok,Qtmp,statusTmp] = fixed_phase_feasibility( ...
        high,Nt,K,P0,Hscaled,noiseScaled,gUsers,ch.sigma_c2,opts, ...
        cvxQuiet,precisionMode);
    numFeasibilityChecks = numFeasibilityChecks + 1;
    lastStatus = statusTmp;
    if ok
        low = high;
        bestQ = Qtmp;
        highFeasible = true;
        if high >= upperCap*(1-1e-12)
            break;
        end
        high = min(2*high,upperCap);
    else
        highFeasible = false;
        break;
    end
end

% If the cap itself is feasible, return the cap as a lower-bound estimate.
% Otherwise refine the feasible/infeasible bracket.
if ~highFeasible
    for ib = 1:maxBisection
        mid = 0.5*(low+high);
        [ok,Qtmp,statusTmp] = fixed_phase_feasibility( ...
            mid,Nt,K,P0,Hscaled,noiseScaled,gUsers,ch.sigma_c2,opts, ...
            cvxQuiet,precisionMode);
        numFeasibilityChecks = numFeasibilityChecks + 1;
        lastStatus = statusTmp;
        if ok
            low = mid;
            bestQ = Qtmp;
        else
            high = mid;
        end
        if (high-low)/max(low,realmin) <= relTol
            break;
        end
    end
end

gammaMax = low;
sinr = fixed_phase_sinr(bestQ,gUsers,ch.sigma_c2,K);
sol = make_solution(true,gammaMax,bestQ,sinr,'fixed-phase-pure-comm-solved');
info = make_info(true,numFeasibilityChecks,lastStatus,low,high,precisionMode);
end

function [ok,Qcomm,statusText] = fixed_phase_feasibility( ...
    gamma,Nt,K,P0,Hscaled,noiseScaled,gUsers,sigmaC2,opts, ...
    cvxQuiet,precisionMode)
        Qcomm = zeros(Nt,Nt,K);
        statusText = 'not-solved';
        try
            if cvxQuiet
                cvx_begin sdp quiet
            else
                cvx_begin sdp
            end
                if isfield(opts,'cvxSolver') && ~isempty(opts.cvxSolver)
                    cvx_solver(char(string(opts.cvxSolver)));
                end
                switch precisionMode
                    case 'default'
                        cvx_precision default
                    case 'high'
                        cvx_precision high
                    case 'best'
                        cvx_precision best
                end
                % Some CVX/MATLAB combinations do not instantiate a 3-D
                % Hermitian-semidefinite variable correctly.  One PSD block
                % matrix is equivalent here because only its diagonal user
                % blocks enter the objective and constraints.
                variable Xbig(Nt*K,Nt*K) hermitian semidefinite
                expression Xtotal(Nt,Nt)
                Xtotal = Xbig(1:Nt,1:Nt);
                for jj = 2:K
                    idxj = (jj-1)*Nt+(1:Nt);
                    Xtotal = Xtotal + Xbig(idxj,idxj);
                end
                minimize(0)
                subject to
                    real(trace(Xtotal)) <= 1;
                    for kk = 1:K
                        idxk = (kk-1)*Nt+(1:Nt);
                        real(trace(Hscaled(:,:,kk)*Xbig(idxk,idxk))) >= ...
                            gamma*(real(trace(Hscaled(:,:,kk)*( ...
                            Xtotal-Xbig(idxk,idxk)))) + noiseScaled(kk));
                    end
            cvx_end
        catch ME
            try
                cvx_clear;
            catch
            end
            ok = false;
            statusText = ['exception: ',ME.message];
            return;
        end
        statusText = char(string(cvx_status));
        ok = contains(lower(statusText),'solved');
        if ~ok
            return;
        end
        for kk = 1:K
            idxk = (kk-1)*Nt+(1:Nt);
            Qcomm(:,:,kk) = P0*full(Xbig(idxk,idxk));
            Qcomm(:,:,kk) = (Qcomm(:,:,kk)+Qcomm(:,:,kk)')/2;
        end
        sinrCheck = fixed_phase_sinr(Qcomm,gUsers,sigmaC2,K);
        tol = 1e-8 + 1e-5*gamma;
        ok = all(sinrCheck + tol >= gamma);
        if ~ok
            statusText = [statusText,'; exact-SINR-postcheck-failed'];
        end
end

function sinr = fixed_phase_sinr(Qcomm,gUsers,sigmaC2,K)
        sinr = zeros(K,1);
        for kk = 1:K
            gk = gUsers(:,kk);
            desired = real(gk'*Qcomm(:,:,kk)*gk);
            interference = 0;
            for jj = 1:K
                if jj ~= kk
                    interference = interference + real(gk'*Qcomm(:,:,jj)*gk);
                end
            end
            sinr(kk) = desired/max(interference+sigmaC2,realmin);
        end
end

function sol = make_solution(success,gammaMax,Qcomm,sinr,status)
sol = struct( ...
    'success',logical(success), ...
    'gammaMax',gammaMax, ...
    'gammaMaxDB',10*log10(max(gammaMax,realmin)), ...
    'Qcomm',Qcomm, ...
    'sinr',sinr, ...
    'status',status);
end

function info = make_info(success,numChecks,lastStatus,low,high,precisionMode)
info = struct( ...
    'success',logical(success), ...
    'numFeasibilityChecks',numChecks, ...
    'lastCVXStatus',lastStatus, ...
    'lowerBound',low, ...
    'upperBound',high, ...
    'precisionMode',precisionMode, ...
    'status','fixed-phase-pure-communication-bound');
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
