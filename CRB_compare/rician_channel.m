function H = rician_channel(rxPos, txPos, params, seed)
%RICIAN_CHANNEL Far-field Rician channel with steering-vector LoS component.
%
% This function is used only for far-field cascaded segments, e.g., the
% outermost SIM layer to CUs/targets and all multi-hop rRIS links. It should
% NOT be used for SIM internal near-field propagation; the SIM stack uses
% diffraction_matrix.m instead.
%
% Model:
%   H = sqrt(PL(d0)) * ( sqrt(kappa/(kappa+1))*H_LoS
%                      + sqrt(1/(kappa+1))*H_NLoS ),
% where d0 is the center-to-center distance and H_NLoS has i.i.d. CN(0,1)
% entries. The LoS part uses far-field array steering vectors, not element-
% wise spherical-wave distances.

if nargin >= 4 && ~isempty(seed)                    % If a seed is provided for this link,
    oldStream = rng;                                % store the current random state.
    rng(seed);                                      % use a deterministic NLoS realization for this link.
else                                                % Otherwise,
    oldStream = [];                                 % do not restore any state later.
end

if size(rxPos,2) ~= 3                               % If receiver coordinate is given as a vector,
    rxPos = reshape(rxPos,1,3);                     % convert it to one 3-D receiver coordinate.
end
if size(txPos,2) ~= 3                               % If transmitter coordinate is given as a vector,
    txPos = reshape(txPos,1,3);                     % convert it to one 3-D transmitter coordinate.
end

Nr = size(rxPos,1);                                 % Number of receive points/elements.
Nt = size(txPos,1);                                 % Number of transmit points/elements.
lambda = params.lambda;                            % Carrier wavelength.
rxCenter = mean(rxPos,1);                          % Receiver array center.
txCenter = mean(txPos,1);                          % Transmitter array center.
linkVec = rxCenter - txCenter;                     % Direction from TX center to RX center.
d0 = norm(linkVec);                                % Center-to-center link distance.
u = linkVec / max(d0,1e-15);                       % Unit propagation direction.
pl = pathloss_linear(d0, params);                  % Large-scale pathloss based on center distance.

rxRel = rxPos - rxCenter;                          % Receiver element positions relative to RX center.
txRel = txPos - txCenter;                          % Transmitter element positions relative to TX center.
a_rx = exp(-1j*2*pi/lambda * (rxRel*u.'));         % RX far-field steering response, Nr-by-1.
a_tx = exp(-1j*2*pi/lambda * (txRel*u.'));         % TX far-field steering response, Nt-by-1.
Hlos = a_rx * a_tx';                               % LoS component a_rx a_tx^H, size Nr-by-Nt.
Hnlos = (randn(Nr,Nt)+1j*randn(Nr,Nt))/sqrt(2);    % NLoS component with CN(0,1) entries.

kappa = params.ricianK;                            % Linear Rician factor, not dB.
Hsmall = sqrt(kappa/(kappa+1))*Hlos ...            % Weighted LoS component.
       + sqrt(1/(kappa+1))*Hnlos;                  % Weighted random NLoS component.
H = sqrt(pl) * Hsmall;                             % Apply large-scale pathloss to the small-scale channel.

if ~isempty(oldStream)                              % If this function changed the random seed,
    rng(oldStream);                                 % restore the previous random state.
end
end
