function bank = build_common_rRIS_bank(params, Lmax)
%BUILD_COMMON_RRIS_BANK Generate fixed rRIS geometry and channel bank once.
%
% Purpose:
%   For a fair comparison of different rRIS hop numbers, all physical links
%   that are common across L should be generated once outside the L-loop and
%   then reused. RIS positions are fixed. The first RIS remains at the same
%   position for all L, and users/targets are placed at the same relative
%   distance and direction with respect to the last active RIS.
%
%   L=1 uses RIS_1 as the last RIS.
%   L=2 uses RIS_1 -> RIS_2, with RIS_2 as the last RIS.
%   L=3 uses RIS_1 -> RIS_2 -> RIS_3, with RIS_3 as the last RIS.
%
% Output fields:
%   bank.bsPos       : fixed BS array coordinates.
%   bank.layerPos    : fixed coordinates of RIS_1,...,RIS_Lmax.
%   bank.G           : fixed BS-to-RIS_1 channel.
%   bank.Omega       : fixed RIS_{ell-1}-to-RIS_ell channels.
%   bank.userOffset  : user positions relative to the last active RIS center.
%   bank.targetOffset: target positions relative to the last active RIS center.
%   bank.hUsersOuter : reusable last-RIS-to-user channels in local geometry.
%   bank.HsOuter/RhOuter : reusable target response in local geometry.

if nargin < 2 || isempty(Lmax)
    Lmax = 3;
end

Nt = params.Nt;                                    % Number of BS antennas.
N = params.N;                                      % Number of elements per rRIS layer.
K = params.K;                                      % Number of communication users.
d = params.d;                                      % Element spacing.

bsPos = upa_positions(Nt, d, [0,0,0], 'yz');       % Fixed BS UPA at x=0.
layerPos = cell(Lmax,1);                           % Fixed rRIS positions.
for ell = 1:Lmax
    x = ell * params.rris.hopDistance;             % RIS_ell center location.
    layerPos{ell} = upa_positions(N, d, [x,0,0], 'yz');
end

G = rician_channel(layerPos{1}, bsPos, params, params.seed.rrisBase + 1);
Omega = cell(Lmax,1);
for ell = 2:Lmax
    Omega{ell} = rician_channel(layerPos{ell}, layerPos{ell-1}, params, params.seed.rrisBase + ell);
end

% Build a local reference last-RIS geometry centered at the origin. Because
% all rRIS layers have the same orientation and aperture, a channel generated
% in this local geometry can be reused for the last active RIS of any L.
refLayerPos = upa_positions(N, d, [0,0,0], 'yz');
refCenter = mean(refLayerPos,1);

userRefPos = position_from_angles(refCenter, ...
    params.rris.lastHopUserDistance, params.rris.userAzDeg, params.rris.userElDeg);
targetRefPos = position_from_angles(refCenter, ...
    params.rris.lastHopTargetDistance, params.rris.targetAzDeg, params.rris.targetElDeg);

userOffset = userRefPos - refCenter;               % Reused relative user geometry.
targetOffset = targetRefPos - refCenter;           % Reused relative target geometry.

hUsersOuter = cell(K,1);                           % Reusable last-RIS-to-CU channels.
for k = 1:K
    Htmp = rician_channel(userRefPos(k,:), refLayerPos, params, params.seed.userBase + k);
    hUsersOuter{k} = Htmp.';                       % Store h_k as N-by-1.
end

paramsTarget = params;                             % Local copy used only to fix target coefficients.
oldStream = rng;                                   % Preserve global random state.
rng(params.seed.targetBase);                       % One target small-scale realization for all L.
paramsTarget.target.xi = (randn(params.Qtar,1) + 1j*randn(params.Qtar,1))/sqrt(2);
rng(oldStream);                                    % Restore global random state.

[HsOuter, RhOuter] = build_target_response(refLayerPos, targetRefPos, paramsTarget);

bank.Lmax = Lmax;                                  % Maximum number of rRIS layers in this bank.
bank.bsPos = bsPos;                                % Fixed BS coordinates.
bank.layerPos = layerPos;                          % Fixed rRIS coordinates.
bank.G = G;                                        % Fixed BS-to-RIS_1 channel.
bank.Omega = Omega;                                % Fixed inter-rRIS channels.
bank.userOffset = userOffset;                      % User positions relative to last active RIS.
bank.targetOffset = targetOffset;                  % Target positions relative to last active RIS.
bank.hUsersOuter = hUsersOuter;                    % Common last-hop user channels.
bank.HsOuter = HsOuter;                            % Common local target response.
bank.RhOuter = RhOuter;                            % Common local target covariance prior.
bank.targetXi = paramsTarget.target.xi;            % Stored target coefficients for reproducibility.
end
