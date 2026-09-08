function out = SIM_solve_comm_margin_SDR(params,ch,theta,gammaTarget,alg)
%SIM_SOLVE_COMM_MARGIN_SDR Covariance-domain communication Phase-I SDP.
% The hidden carrier uses lifted covariances only; no rank-one recovery.

K=params.K; Nt=params.Nt; P0=params.P0;
gammaTarget=gammaTarget(:);
% Solve the hidden carrier against a slightly stricter target, then assess
% it against the requested target.  This absorbs CVX boundary/backward
% error without relaxing any formal SINR constraint.
guardRel=max(0,get_option(alg,'carrierSDRSINRGuardRel',5e-4));
gammaConstraint=gammaTarget*(1+guardRel);
[P,~]=build_P(theta,ch.Omega);
template=struct('success',false,'status','not-run','cvxStatus','not-run', ...
    'Qcomm',zeros(Nt,Nt,K),'Rs',zeros(Nt),'Rx',zeros(Nt), ...
    'sinr',NaN(K,1),'minRelativeMargin',-Inf,'scaledResidual',NaN(K,1), ...
    'covarianceFloorRel',NaN,'power',NaN);
out=template;

Hscaled=zeros(Nt,Nt,K); noiseScaled=zeros(K,1);
for k=1:K
    gk=ch.G'*P'*ch.hUsers{k}; Hk=gk*gk'; Hk=(Hk+Hk')/2;
    userScale=max(P0*real(trace(Hk))+ch.sigma_c2,realmin);
    Hscaled(:,:,k)=(P0/userScale)*Hk;
    noiseScaled(k)=ch.sigma_c2/userScale;
end

floorList=get_option(alg,'carrierCovarianceFloorRelList',[1e-3,1e-5,0]);
floorList=floorList(:).';
floorList=floorList(isfinite(floorList)&floorList>=0&floorList<1);
if isempty(floorList), floorList=0; end
best=[];
for ii=1:numel(floorList)
    trial=solve_one_floor(params,ch,P,gammaConstraint,gammaTarget, ...
        Hscaled,noiseScaled,floorList(ii),template);
    if trial.success
        if isempty(best)||trial.minRelativeMargin>best.minRelativeMargin
            best=trial;
        end
        if trial.minRelativeMargin>=-get_option(alg,'carrierSINRRelativeTol',2e-5)
            break;
        end
    end
end
if ~isempty(best), out=best; end
end

function trial=solve_one_floor(params,ch,P,gammaConstraint,gammaTarget, ...
    Hscaled,noiseScaled,floorRel,template)
% CVX must run in a non-nested workspace for 3-D Hermitian variables.
K=params.K; Nt=params.Nt; P0=params.P0;
trial=template; trial.covarianceFloorRel=floorRel;
try
    cvx_clear;
catch
end
try
    cvx_begin sdp quiet
        cvx_precision high
        variable Xcomm(Nt,Nt,K) hermitian semidefinite
        variable Xs(Nt,Nt) hermitian semidefinite
        variable t
        expression Xtotal(Nt,Nt)
        Xtotal=Xs;
        for jj=1:K, Xtotal=Xtotal+Xcomm(:,:,jj); end
        maximize(t)
        subject to
            real(trace(Xtotal))<=1;
            Xtotal-floorRel/Nt*eye(Nt)>=0;
            for kk=1:K
                real(trace(Hscaled(:,:,kk)*Xcomm(:,:,kk))) - ...
                    gammaConstraint(kk)*(real(trace(Hscaled(:,:,kk)*( ...
                    Xtotal-Xcomm(:,:,kk))))+noiseScaled(kk))>=t;
            end
    cvx_end
catch ME
    trial.status=['CVX exception: ',ME.message]; trial.cvxStatus='exception';
    return;
end
trial.cvxStatus=char(string(cvx_status));
if ~contains(lower(trial.cvxStatus),'solved')|| ...
        any(~isfinite(Xcomm(:)))||any(~isfinite(Xs(:)))
    trial.status=['margin SDP failed: ',trial.cvxStatus]; return;
end
for kk=1:K
    trial.Qcomm(:,:,kk)=P0*full(Xcomm(:,:,kk));
    trial.Qcomm(:,:,kk)=(trial.Qcomm(:,:,kk)+trial.Qcomm(:,:,kk)')/2;
end
trial.Rs=P0*full(Xs); trial.Rs=(trial.Rs+trial.Rs')/2;
trial.Rx=trial.Rs;
for kk=1:K, trial.Rx=trial.Rx+trial.Qcomm(:,:,kk); end
trial.Rx=(trial.Rx+trial.Rx')/2; trial.power=real(trace(trial.Rx));
trial.sinr=covariance_sinr(params,ch,P,trial.Qcomm,trial.Rs);
trial.minRelativeMargin=min(trial.sinr./max(gammaTarget,realmin)-1);
trial.scaledResidual=zeros(K,1);
for kk=1:K
    XtotalValue=trial.Rx/P0;
    desiredValue=real(trace(Hscaled(:,:,kk)*(trial.Qcomm(:,:,kk)/P0)));
    interferenceValue=real(trace(Hscaled(:,:,kk)*( ...
        XtotalValue-trial.Qcomm(:,:,kk)/P0)))+noiseScaled(kk);
    trial.scaledResidual(kk)=desiredValue-gammaTarget(kk)*interferenceValue;
end
trial.success=all(isfinite(trial.sinr));
trial.status='communication-margin-SDR-solved';
end

function sinr=covariance_sinr(params,ch,P,Qcomm,Rs)
Rx=Rs; for jj=1:params.K, Rx=Rx+Qcomm(:,:,jj); end
sinr=zeros(params.K,1);
for kk=1:params.K
    gk=ch.G'*P'*ch.hUsers{kk};
    desired=real(gk'*Qcomm(:,:,kk)*gk);
    denominator=real(gk'*(Rx-Qcomm(:,:,kk))*gk)+ch.sigma_c2;
    sinr(kk)=desired/max(denominator,realmin);
end
end

function value=get_option(s,name,defaultValue)
value=defaultValue;
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), value=s.(name); end
end
