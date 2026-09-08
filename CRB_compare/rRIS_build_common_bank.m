function bank = rRIS_build_common_bank(params,Lmax)
%RRIS_BUILD_COMMON_BANK Fixed multi-hop rRIS geometry/channel realization.
% Physical construction matches MIc-Mis_revised/build_common_rRIS_bank.m.

if nargin<2 || isempty(Lmax), Lmax=3; end
Nt=params.Nt; N=params.N; K=params.K; d=params.d;
bsPos=upa_positions(Nt,d,[0,0,0],'yz');
layerPos=cell(Lmax,1);
for ell=1:Lmax
    layerPos{ell}=upa_positions(N,d, ...
        [ell*params.rris.hopDistance,0,0],'yz');
end
G=rician_channel(layerPos{1},bsPos,params,params.seed.rrisBase+1);
Omega=cell(Lmax,1);
for ell=2:Lmax
    Omega{ell}=rician_channel(layerPos{ell},layerPos{ell-1}, ...
        params,params.seed.rrisBase+ell);
end

refLayerPos=upa_positions(N,d,[0,0,0],'yz');
refCenter=mean(refLayerPos,1);
userRefPos=position_from_angles(refCenter, ...
    params.rris.lastHopUserDistance,params.rris.userAzDeg, ...
    params.rris.userElDeg);
targetRefPos=position_from_angles(refCenter, ...
    params.rris.lastHopTargetDistance,params.rris.targetAzDeg, ...
    params.rris.targetElDeg);
userOffset=userRefPos-refCenter;
targetOffset=targetRefPos-refCenter;
hUsersOuter=cell(K,1);
for k=1:K
    tmp=rician_channel(userRefPos(k,:),refLayerPos,params, ...
        params.seed.userBase+k);
    hUsersOuter{k}=tmp.';
end

paramsTarget=params;
oldStream=rng;
rng(params.seed.targetBase);
paramsTarget.target.xi=(randn(params.Qtar,1)+ ...
    1j*randn(params.Qtar,1))/sqrt(2);
rng(oldStream);
[HsOuter,RhOuter]=build_target_response( ...
    refLayerPos,targetRefPos,paramsTarget);

bank=struct('Lmax',Lmax,'bsPos',bsPos,'layerPos',{layerPos}, ...
    'G',G,'Omega',{Omega},'userOffset',userOffset, ...
    'targetOffset',targetOffset,'hUsersOuter',{hUsersOuter}, ...
    'HsOuter',HsOuter,'RhOuter',RhOuter, ...
    'targetXi',paramsTarget.target.xi);
end
