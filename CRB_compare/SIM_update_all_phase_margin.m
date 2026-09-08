function [theta,info] = SIM_update_all_phase_margin( ...
    params,ch,theta,Qcomm,Rs,gammaTarget,alg)
%SIM_UPDATE_ALL_PHASE_MARGIN Communication-only ascent over all SIM layers.
%
% This Phase-I helper is allowed to move Phi_1,...,Phi_L because its state is
% a hidden communication carrier, not the reported sensing solution.

maxIter = max(0,round(get_option(alg,'carrierAllPhaseMaxIter',20)));
maxLineSearch = max(1,round(get_option(alg,'carrierPhaseMaxLineSearch',12)));
alpha0 = get_option(alg,'carrierPhaseAlpha0',0.15);
alphaMin = get_option(alg,'carrierPhaseAlphaMin',1e-6);
beta = get_option(alg,'carrierPhaseLineSearchBeta',0.5);
armijoC = get_option(alg,'carrierMarginArmijoC',1e-4);
improveTol = get_option(alg,'carrierMarginImproveTol',1e-10);
activeTol = get_option(alg,'carrierMarginActiveTol',1e-7);
feasTol = get_option(alg,'carrierSINRRelativeTol',2e-5);
requireCRBHealth = logical(get_option(alg,'carrierRequireCRBHealth',true));
maxEigDrop = max(1,get_option(alg,'carrierMaxCRBEigDropFactor',10));

[margin,G] = relative_margin_gradient(theta,true);
minMargin = min(margin);
[healthStart,healthFloorA,healthFloorB] = crb_health(theta,maxEigDrop);
info = struct('success',all(isfinite(margin)),'status','not-started', ...
    'iters',0,'acceptedSteps',0,'marginStart',minMargin, ...
    'marginEnd',minMargin,'lineSearchTrials',0,'lastAlpha',NaN, ...
    'crbHealthStart',healthStart,'crbHealthRejects',0);
if ~info.success || maxIter == 0
    info.status = 'invalid-or-disabled';
    return;
end
if minMargin >= -feasTol
    info.status = 'communication-feasible';
    return;
end

for it = 1:maxIter
    info.iters = it;
    [direction,prediction] = min_margin_direction(margin,G,activeTol);
    if ~isfinite(prediction) || prediction <= 0 || norm(direction) == 0
        info.status = 'no-common-first-order-ascent';
        break;
    end

    accepted = false;
    alpha = alpha0;
    for ls = 1:maxLineSearch
        info.lineSearchTrials = info.lineSearchTrials+1;
        thetaTry = mod(theta+alpha*reshape(direction,size(theta)),2*pi);
        marginTry = relative_margin_gradient(thetaTry,false);
        minTry = min(marginTry);
        healthOK = true;
        if requireCRBHealth
            healthOK = crb_health(thetaTry,healthFloorA,healthFloorB);
            if ~healthOK, info.crbHealthRejects=info.crbHealthRejects+1; end
        end
        if healthOK && all(isfinite(marginTry)) && ...
                minTry >= minMargin+armijoC*alpha*prediction-improveTol
            accepted = true;
            break;
        end
        alpha = beta*alpha;
        if alpha < alphaMin, break; end
    end
    if ~accepted
        info.status = 'line-search-failed';
        break;
    end

    theta = thetaTry;
    margin = marginTry;
    minMargin = minTry;
    info.acceptedSteps = info.acceptedSteps+1;
    info.lastAlpha = alpha;
    if minMargin >= -feasTol
        info.status = 'communication-feasible';
        break;
    end
    [margin,G] = relative_margin_gradient(theta,true);
    minMargin = min(margin);
end
info.marginEnd = minMargin;
info.success = isfinite(minMargin);
if strcmp(info.status,'not-started')
    info.status = 'maximum-iterations-reached';
end

    function [relativeMargin,Gtheta] = relative_margin_gradient(thetaValue,needGradient)
        [P,~] = build_P(thetaValue,ch.Omega);
        Rx = Rs;
        for jj = 1:params.K, Rx=Rx+Qcomm(:,:,jj); end
        relativeMargin = zeros(params.K,1);
        if needGradient
            Gtheta = zeros(numel(thetaValue),params.K);
        else
            Gtheta = [];
        end
        for kk = 1:params.K
            gk = ch.G'*P'*ch.hUsers{kk};
            desired = real(gk'*Qcomm(:,:,kk)*gk);
            denominator = real(gk'*(Rx-Qcomm(:,:,kk))*gk)+ch.sigma_c2;
            denominator = max(denominator,realmin);
            relativeMargin(kk) = desired/(max(gammaTarget(kk),realmin)*denominator)-1;
            if ~needGradient, continue; end
            for ell = 1:ch.L
                [UL,UR] = compute_UL_UR(thetaValue,ch.Omega,ell);
                for nn = 1:params.N
                    En = zeros(params.N); En(nn,nn)=1;
                    dP = 1j*exp(1j*thetaValue(nn,ell))*UL*En*UR;
                    dg = ch.G'*dP'*ch.hUsers{kk};
                    dDesired = 2*real(dg'*Qcomm(:,:,kk)*gk);
                    dTotal = 2*real(dg'*Rx*gk);
                    dDenominator = dTotal-dDesired;
                    dMargin = (dDesired*denominator-desired*dDenominator) / ...
                        (max(gammaTarget(kk),realmin)*denominator^2);
                    Gtheta(nn+(ell-1)*params.N,kk) = real(dMargin);
                end
            end
        end
    end

    function [ok,floorA,floorB] = crb_health(thetaValue,arg2,arg3)
        [Pvalue,~] = build_P(thetaValue,ch.Omega);
        RxValue = get_option(alg,'carrierCRBHealthRx',[]);
        if isempty(RxValue)
            RxValue = Rs;
            for jj=1:params.K, RxValue=RxValue+Qcomm(:,:,jj); end
        end
        [~,crbInfo] = SIM_crb_value(ch,Pvalue,RxValue,alg);
        relA = get_option(crbInfo,'relMinEigA',-Inf);
        relB = get_option(crbInfo,'relMinEigB',-Inf);
        if nargin==2
            dropFactor=arg2;
            baseFloor=max(get_option(alg,'carrierMinCRBRelativeEig',1e-12),0);
            floorA=max(baseFloor,relA/dropFactor);
            floorB=max(baseFloor,relB/dropFactor);
        else
            floorA=arg2; floorB=arg3;
        end
        ok=logical(get_option(crbInfo,'valid',false)) && ...
            isfinite(relA) && isfinite(relB) && relA>=floorA && relB>=floorB;
    end
end

function [direction,prediction] = min_margin_direction(margin,G,activeTol)
minimum = min(margin);
active = find(margin <= minimum+activeTol);
if isempty(active), [~,active]=min(margin); end
if numel(active)==1
    gbar = G(:,active);
elseif numel(active)==2
    g1=G(:,active(1)); g2=G(:,active(2)); v=g1-g2;
    vv=real(v'*v);
    if vv<=realmin, lambda=0.5;
    else, lambda=min(max(real(g2'*(g2-g1))/vv,0),1); end
    gbar=lambda*g1+(1-lambda)*g2;
else
    gbar=mean(G(:,active),2);
end
ng=norm(gbar);
if ~isfinite(ng) || ng==0
    direction=zeros(size(gbar)); prediction=0; return;
end
direction=gbar/ng;
prediction=min(real(G(:,active)'*direction));
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
