function metrics = SIM_evaluate_CRB_SINR_metrics( ...
    params, ch, P, Qcomm, Rs, S, Rx, gammaTarget, alg) %#ok<INUSD>
%SIM_EVALUATE_CRB_SINR_METRICS Evaluate an SDR covariance solution.
% No beamformer/rank-one recovery is performed here.
K = params.K; gammaTarget = gammaTarget(:);
sinr=zeros(K,1); rawConstraint=zeros(K,1); scaledConstraint=zeros(K,1);
MIc=0; rate=0;
for k=1:K
    gk=ch.G'*P'*ch.hUsers{k};
    desired=real(gk'*Qcomm(:,:,k)*gk);
    interfNoise=max(real(gk'*Rx*gk)+ch.sigma_c2-desired,1e-30);
    sinr(k)=desired/interfNoise;
    rawConstraint(k)=gammaTarget(k)*interfNoise-desired;
    scale=max(abs(gammaTarget(k)*interfNoise)+abs(desired),1e-30);
    scaledConstraint(k)=rawConstraint(k)/scale;
    MIc=MIc+log(1+sinr(k)); rate=rate+log2(1+sinr(k));
end
[Jcrb,crbInfo]=SIM_crb_value(ch,P,Rx,alg);
relTol=get_option(alg,'sdrPostSINRRelTol',1e-5);
absTol=get_option(alg,'sdrPostSINRAbsTol',1e-10);
sinrTolerance=absTol+relTol.*abs(gammaTarget);
metrics.CRB=Jcrb; metrics.crbInfo=crbInfo;
metrics.traceAinv=crbInfo.traceAinv; metrics.traceBinv=crbInfo.traceBinv;
metrics.logCRB=crbInfo.logCRB; metrics.relMinEigA=crbInfo.relMinEigA;
metrics.relMinEigB=crbInfo.relMinEigB; metrics.conditionWarning=crbInfo.conditionWarning;
metrics.lambdaMinA=crbInfo.minEigA; metrics.rcondA=crbInfo.rcondA;
metrics.rankA=NaN; metrics.condA=1/max(crbInfo.rcondA,realmin);
metrics.lambdaMinB=crbInfo.minEigB; metrics.rcondB=crbInfo.rcondB;
metrics.rankB=NaN; metrics.condB=1/max(crbInfo.rcondB,realmin);
metrics.MIc=MIc; metrics.rate=rate;
% Waveform-dependent diagnostics are deferred because forming them would
% require selecting a factor of Rx during the optimization.
metrics.MIs=NaN; metrics.nmmse=NaN;
metrics.sinr=sinr; metrics.rawConstraint=rawConstraint;
metrics.scaledConstraint=scaledConstraint;
metrics.maxScaledConstraintViolation=max([0;scaledConstraint]);
metrics.power=real(trace(Rx)); metrics.sinrTolerance=sinrTolerance;
metrics.sinrFeasible=all(sinr+sinrTolerance>=gammaTarget);
metrics.crbValid=crbInfo.valid; metrics.Qcomm=Qcomm; metrics.Rs=Rs;
end

function value=get_option(s,name,defaultValue)
value=defaultValue;
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)), value=s.(name); end
end
