function common = build_common_SIM_links(params)
%BUILD_COMMON_SIM_LINKS Generate SIM outermost-layer user/target links once.
%
% Purpose:
%   In the SIM comparison, the BS-to-outermost-layer distance is fixed as
%   params.sim.totalThickness for all L. Therefore, the outermost-layer-to-
%   user channels and target response should be generated once outside the
%   L-loop and reused for L=1,2,3. This avoids changing the random channel
%   realization when only the number of SIM layers is changed.
%
% Output fields:
%   common.outerLayerPos : coordinates of the common outermost SIM layer.
%   common.userPos       : user positions relative to the common outer layer.
%   common.targetPos     : target positions relative to the common outer layer.
%   common.hUsers        : fixed outermost-layer-to-user channels.
%   common.Hs, common.Rh : fixed target response and covariance prior.

N = params.N;                                      % Number of SIM elements per layer.
K = params.K;                                      % Number of communication users.
d = params.d;                                      % Element spacing.

outerLayerPos = upa_positions(N, d, ...            % Common outermost SIM layer for every L.
    [params.sim.totalThickness, 0, 0], 'yz');
outerCenter = mean(outerLayerPos,1);               % Center of the outermost SIM layer.

userPos = position_from_angles(outerCenter, ...    % Fixed user positions relative to outer layer.
    params.sim.userDistance, params.sim.userAzDeg, params.sim.userElDeg);
targetPos = position_from_angles(outerCenter, ...  % Fixed target positions relative to outer layer.
    params.sim.targetDistance, params.sim.targetAzDeg, params.sim.targetElDeg);

hUsers = cell(K,1);                                % Fixed outermost-layer-to-CU channels.
for k = 1:K
    Htmp = rician_channel(userPos(k,:), outerLayerPos, params, params.seed.userBase + k);
    hUsers{k} = Htmp.';                            % Store h_k as N-by-1.
end

paramsTarget = params;                             % Local copy used only to fix target coefficients.
oldStream = rng;                                   % Preserve global random state.
rng(params.seed.targetBase);                       % One target small-scale realization for this scenario.
paramsTarget.target.xi = (randn(params.Qtar,1) + 1j*randn(params.Qtar,1))/sqrt(2);
rng(oldStream);                                    % Restore global random state.

[Hs, Rh] = build_target_response(outerLayerPos, targetPos, paramsTarget);

common.outerLayerPos = outerLayerPos;              % Common outermost layer coordinates.
common.userPos = userPos;                          % Common user positions.
common.targetPos = targetPos;                      % Common target positions.
common.hUsers = hUsers;                            % Common user channels.
common.Hs = Hs;                                    % Common instantaneous target response.
common.Rh = Rh;                                    % Common target covariance prior.
common.targetXi = paramsTarget.target.xi;          % Stored target coefficients for reproducibility.
end
