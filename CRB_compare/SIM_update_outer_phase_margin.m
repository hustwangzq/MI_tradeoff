function [theta,info] = SIM_update_outer_phase_margin( ...
    params,ch,theta,Qcomm,Rs,gammaTarget,alg)
%SIM_UPDATE_OUTER_PHASE_MARGIN Improve communication reserve with Phi_L.
%
% For the extended-target CRB used by this project, the outermost diagonal
% phase Phi_L acts by unitary congruence on both sensing FIM factors. Hence
% it is CRB-invariant and is optimized here only for the true covariance-
% domain SINR margin. No rank-one recovery is performed.

gammaTarget = gammaTarget(:);
L = ch.L;
maxIter = max(0,round(get_option(alg,'outerMarginMaxIter',12)));
maxLineSearch = max(1,round(get_option(alg,'outerMarginMaxLineSearch',10)));
alpha0 = get_option(alg,'outerMarginAlpha0',0.20);
alphaMin = get_option(alg,'outerMarginAlphaMin',1e-5);
beta = get_option(alg,'outerMarginLineSearchBeta',0.5);
smoothTau = get_option(alg,'outerMarginSmoothTau',100);
gradTol = get_option(alg,'outerMarginGradTol',1e-8);
improveTol = get_option(alg,'outerMarginImproveTol',1e-7);
crbRelTol = get_option(alg,'outerMarginCRBRelTol',1e-7);

Rx = Rs;
for kk = 1:params.K
    Rx = Rx + Qcomm(:,:,kk);
end
Rx = (Rx+Rx')/2;

[P0,~] = build_P(theta,ch.Omega);
[crb0,crbInfo0] = SIM_crb_value(ch,P0,Rx,alg);
[margin0,~,sinr0] = margin_value_gradient( ...
    params,ch,theta,Qcomm,Rs,gammaTarget,L,smoothTau,false);

info.success = crbInfo0.valid && isfinite(margin0);
info.status = 'not-started';
info.iters = 0;
info.acceptedSteps = 0;
info.marginStart = margin0;
info.marginEnd = margin0;
info.minSINRStart = min(sinr0);
info.minSINREnd = info.minSINRStart;
info.lastAlpha = NaN;
info.lastGradNorm = NaN;
info.crbRelativeDriftMax = 0;

if ~logical(get_option(alg,'outerMarginEnable',true)) || maxIter == 0
    info.status = 'outer-margin update disabled';
    return;
end
if ~info.success
    info.status = 'invalid initial outer-margin state';
    return;
end

for it = 1:maxIter
    info.iters = it;
    [marginCur,gradOuter,sinrCur] = margin_value_gradient( ...
        params,ch,theta,Qcomm,Rs,gammaTarget,L,smoothTau,true);
    gradNorm = norm(gradOuter);
    info.lastGradNorm = gradNorm;
    if ~isfinite(gradNorm) || gradNorm <= gradTol
        info.status = 'outer-margin gradient tolerance reached';
        break;
    end

    direction = gradOuter/max(gradNorm,eps);
    alpha = alpha0;
    accepted = false;
    bestTheta = theta;
    bestMargin = marginCur;
    bestSinr = sinrCur;
    bestDrift = 0;

    for ls = 1:maxLineSearch
        if alpha < alphaMin
            break;
        end
        thetaTry = theta;
        thetaTry(:,L) = mod(theta(:,L)+alpha*direction,2*pi);
        [marginTry,~,sinrTry] = margin_value_gradient( ...
            params,ch,thetaTry,Qcomm,Rs,gammaTarget,L,smoothTau,false);
        [Ptry,~] = build_P(thetaTry,ch.Omega);
        [crbTry,crbInfoTry] = SIM_crb_value(ch,Ptry,Rx,alg);
        drift = abs(crbTry-crb0)/max(abs(crb0),realmin);
        if crbInfoTry.valid && isfinite(marginTry) && ...
                marginTry >= marginCur+improveTol && drift <= crbRelTol
            accepted = true;
            bestTheta = thetaTry;
            bestMargin = marginTry;
            bestSinr = sinrTry;
            bestDrift = drift;
            break;
        end
        alpha = beta*alpha;
    end

    if ~accepted
        info.status = 'outer-margin line search failed';
        break;
    end

    theta = bestTheta;
    info.acceptedSteps = info.acceptedSteps+1;
    info.lastAlpha = alpha;
    info.marginEnd = bestMargin;
    info.minSINREnd = min(bestSinr);
    info.crbRelativeDriftMax = max(info.crbRelativeDriftMax,bestDrift);
    if bestMargin-marginCur <= improveTol
        info.status = 'outer-margin improvement tolerance reached';
        break;
    end
end

if strcmp(info.status,'not-started')
    info.status = 'outer-margin maximum iterations reached';
end
end

function [minMargin,gradOuter,sinr] = margin_value_gradient( ...
    params,ch,theta,Qcomm,Rs,gammaTarget,outerLayer,smoothTau,needGradient)
% Return the exact minimum relative SINR margin. The gradient is that of a
% smooth-min surrogate and is used only to generate an ascent direction.

[P,~] = build_P(theta,ch.Omega);
Rx = Rs;
for kk = 1:params.K
    Rx = Rx + Qcomm(:,:,kk);
end
Rx = (Rx+Rx')/2;

desired = zeros(params.K,1);
denom = zeros(params.K,1);
sinr = zeros(params.K,1);
for k = 1:params.K
    gk = ch.G'*P'*ch.hUsers{k};
    desired(k) = real(gk'*Qcomm(:,:,k)*gk);
    total = real(gk'*Rx*gk)+ch.sigma_c2;
    denom(k) = max(total-desired(k),1e-30);
    sinr(k) = desired(k)/denom(k);
end
relativeMargin = sinr./max(gammaTarget,realmin)-1;
minMargin = min(relativeMargin);
gradOuter = zeros(params.N,1);
if ~needGradient
    return;
end

shift = min(relativeMargin);
weights = exp(-smoothTau*(relativeMargin-shift));
weights = weights/max(sum(weights),realmin);
[UL,UR] = compute_UL_UR(theta,ch.Omega,outerLayer);

for n = 1:params.N
    En = zeros(params.N,params.N);
    En(n,n) = 1;
    dP = 1j*exp(1j*theta(n,outerLayer))*UL*En*UR;
    gradValue = 0;
    for k = 1:params.K
        hk = ch.hUsers{k};
        gk = ch.G'*P'*hk;
        dg = ch.G'*dP'*hk;
        dDesired = 2*real(dg'*Qcomm(:,:,k)*gk);
        dTotal = 2*real(dg'*Rx*gk);
        dDenom = dTotal-dDesired;
        dMargin = (dDesired*denom(k)-desired(k)*dDenom) / ...
            (max(gammaTarget(k),realmin)*denom(k)^2);
        gradValue = gradValue+weights(k)*dMargin;
    end
    gradOuter(n) = real(gradValue);
end
end

function value = get_option(s,name,defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
end
end
