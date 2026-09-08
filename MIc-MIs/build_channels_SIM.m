function ch = build_channels_SIM(params, L, mode, common)
%BUILD_CHANNELS_SIM Build the transmitter-side SIM cascaded channel.
% There is no direct BS-user or BS-target link. The communication path is
% BS -> SIM stack -> CU, and the sensing echo path is
% BS -> SIM stack -> target -> SIM stack -> BS.
%
% Propagation rule:
%   1) BS-to-first-SIM-layer and SIM inter-layer links use the near-field
%      Rayleigh-Sommerfeld diffraction matrix.
%   2) Outermost-SIM-layer-to-CU/target links use fixed far-field Rician
%      channels generated outside the L-loop when common is provided.
%
% Fair-comparison rule:
%   The SIM total thickness is fixed, so the outermost layer has the same
%   position for all L. Therefore, common.hUsers, common.Hs, and common.Rh
%   are reused for L=1,2,3 to avoid changing the random realization when
%   only the number of SIM layers is changed.

if nargin < 4 || isempty(common)                   % Backward compatibility if no common links are supplied.
    common = build_common_SIM_links(params);        % Generate the common outer links once locally.
end

Nt = params.Nt;                                    % Number of BS antennas.
N = params.N;                                      % Number of elements per SIM layer.
lambda = params.lambda;                            % Wavelength.
d = params.d;                                      % Element spacing.

bsPos = upa_positions(Nt, d, [0,0,0], 'yz');        % BS UPA on the yz plane at x=0.
layerX = (1:L) * params.sim.totalThickness / L;     % Fixed outermost distance: x_L = 3 lambda for all L.
layerPos = cell(L,1);                               % Coordinates of all SIM layers.
for ell = 1:L                                      % Place layers inside the fixed SIM thickness.
    layerPos{ell} = upa_positions(N, d, [layerX(ell),0,0], 'yz');
end

if mode ~= 1                                       % SIM stack must be modeled as near-field diffraction.
    error('SIM mode should be 1: Rayleigh-Sommerfeld near-field propagation.');
end
G = diffraction_matrix(layerPos{1}, bsPos, lambda, d^2); % G_st: BS-to-first-layer near-field matrix.

Omega = cell(L,1);                                 % Inter-layer propagation matrices Omega_ell.
for ell = 2:L                                      % For ell>=2, Omega{ell} maps layer ell-1 to layer ell.
    Omega{ell} = diffraction_matrix(layerPos{ell}, layerPos{ell-1}, lambda, d^2); % Near-field inter-layer matrix.
end

% Reuse the common outermost-layer links and target prior. For SIM, the
% outermost layer is at the same position for all L, so the absolute user and
% target positions are also fixed across the L comparison.
hUsers = common.hUsers;                            % Fixed outermost-SIM-to-user channels.
Hs = common.Hs;                                    % Fixed instantaneous target response.
Rh = common.Rh;                                    % Fixed target covariance prior.
userPos = common.userPos;                          % Fixed user positions.
targetPos = common.targetPos;                      % Fixed target positions.

ch.type = 'SIM';                                    % Scenario label.
ch.L = L;                                           % Number of SIM layers.
ch.mode = mode;                                     % Propagation mode.
ch.G = G;                                           % G_st: BS-to-first-layer channel, size N-by-Nt.
ch.Omega = Omega;                                   % Inter-layer propagation matrices.
ch.hUsers = hUsers;                                 % User channels h_k from outermost layer.
ch.Hs = Hs;                                         % Target response matrix H_s.
ch.Rh = Rh;                                         % Prior covariance R_{h_s}.
ch.bsPos = bsPos;                                   % BS element positions.
ch.layerPos = layerPos;                             % SIM layer element positions.
ch.userPos = userPos;                               % User positions.
ch.targetPos = targetPos;                           % Target positions.
ch.totalThickness = params.sim.totalThickness;      % BS-to-outermost-layer distance.
ch.commonOuterLayerPos = common.outerLayerPos;      % Common layer used to generate outer links.
ch.rRayleighLayer = rayleigh_distance(layerPos{L}, lambda); % Far-field boundary of the outermost layer.
end
