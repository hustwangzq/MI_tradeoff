function ch = build_channels_rRIS(params, L, bank)
%BUILD_CHANNELS_RRIS Build the multi-hop reflective RIS cascaded channel.
% There is no direct BS-user or BS-target link. The communication path is
% BS -> rRIS_1 -> ... -> rRIS_L -> CU, and the sensing echo follows the same
% rRIS cascade before and after target reflection.
%
% Fair-comparison rule:
%   The rRIS positions and common physical links are generated once in
%   build_common_rRIS_bank.m and reused for all L. The first rRIS stays at
%   the same position for all L. Users and targets are placed at the same
%   relative distance and direction with respect to the last active rRIS.
%   Therefore, the last-rRIS-to-user channels and target prior are reused in
%   local geometry instead of being randomly regenerated for each L.

if nargin < 3 || isempty(bank)                     % Backward compatibility if no common bank is supplied.
    bank = build_common_rRIS_bank(params, L);       % Generate a local bank up to this L.
end
if L > bank.Lmax                                   % Safety check.
    error('Requested L=%d exceeds the common rRIS bank size Lmax=%d.', L, bank.Lmax);
end

lambda = params.lambda;                            % Carrier wavelength.

bsPos = bank.bsPos;                                % Fixed BS UPA coordinates.
layerPos = bank.layerPos(1:L);                     % Active rRIS layers for this L.

G = bank.G;                                        % Fixed BS-to-rRIS_1 channel.
Omega = cell(L,1);                                 % Active inter-rRIS matrices.
for ell = 2:L
    Omega{ell} = bank.Omega{ell};                  % Reuse fixed rRIS_{ell-1}->rRIS_ell channel.
end

lastCenter = mean(layerPos{L},1);                  % Center of the last active rRIS.
userPos = bsxfun(@plus, bank.userOffset, lastCenter);   % Users keep same relative geometry to the last rRIS.
targetPos = bsxfun(@plus, bank.targetOffset, lastCenter); % Targets keep same relative geometry to the last rRIS.

hUsers = bank.hUsersOuter;                         % Reuse local last-rRIS-to-user channels.
Hs = bank.HsOuter;                                 % Reuse local target response.
Rh = bank.RhOuter;                                 % Reuse local target covariance prior.

ch.type = 'rRIS';                                   % Scenario label.
ch.L = L;                                           % Number of rRIS hops/layers.
ch.mode = 2;                                        % Far-field Rician rRIS mode.
ch.G = G;                                           % G_st: BS-to-first-rRIS channel, size N-by-Nt.
ch.Omega = Omega;                                   % Inter-rRIS matrices.
ch.hUsers = hUsers;                                 % User channels h_k from the last active rRIS.
ch.Hs = Hs;                                         % Target response matrix H_s.
ch.Rh = Rh;                                         % Prior covariance R_{h_s}.
ch.bsPos = bsPos;                                   % BS element positions.
ch.layerPos = layerPos;                             % Active rRIS element positions.
ch.userPos = userPos;                               % User positions for this L.
ch.targetPos = targetPos;                           % Target positions for this L.
ch.coverageDistance = norm(userPos(1,:) - mean(bsPos,1)); % Approximate BS-to-user coverage distance.
ch.rRayleighLayer = rayleigh_distance(layerPos{1}, lambda); % Far-field boundary of each rRIS layer.
ch.bankLmax = bank.Lmax;                            % Size of the common bank used.
end
