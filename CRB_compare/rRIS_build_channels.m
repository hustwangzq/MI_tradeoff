function ch = rRIS_build_channels(params,L,bank)
%RRIS_BUILD_CHANNELS Multi-hop reflective RIS cascaded channel.
% There is no direct BS-user or BS-target path. The common bank is reused
% across L=1,2,3 exactly as in the MIc-MIs rRIS experiment.

if nargin<3 || isempty(bank), bank=rRIS_build_common_bank(params,L); end
assert(L<=bank.Lmax,'Requested L exceeds rRIS bank Lmax.');
layerPos=bank.layerPos(1:L);
Omega=cell(L,1);
for ell=2:L, Omega{ell}=bank.Omega{ell}; end
lastCenter=mean(layerPos{L},1);
userPos=bsxfun(@plus,bank.userOffset,lastCenter);
targetPos=bsxfun(@plus,bank.targetOffset,lastCenter);

ch=struct();
ch.type='rRIS'; ch.L=L; ch.mode=2;
ch.G=bank.G; ch.Omega=Omega;
ch.hUsers=bank.hUsersOuter;
ch.Hs=bank.HsOuter; ch.Rh=bank.RhOuter;
ch.bsPos=bank.bsPos; ch.layerPos=layerPos;
ch.userPos=userPos; ch.targetPos=targetPos;
ch.coverageDistance=norm(userPos(1,:)-mean(bank.bsPos,1));
ch.rRayleighLayer=rayleigh_distance(layerPos{1},params.lambda);
ch.bankLmax=bank.Lmax;
end
