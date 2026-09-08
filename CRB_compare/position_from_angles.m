function pos = position_from_angles(center, distance, azDeg, elDeg)
%POSITION_FROM_ANGLES Generate 3-D points from distance, azimuth, and elevation.
% The broadside direction of the RIS/SIM plane is the positive x-axis. The
% y-axis is the azimuth cross-range direction and the z-axis is elevation.
% Therefore, the actual Euclidean distance from center to each generated
% point is exactly the value specified in distance.

az = azDeg(:) * pi/180;                            % Convert azimuth angles to radians.
el = elDeg(:) * pi/180;                            % Convert elevation angles to radians.
if isscalar(distance)                               % If one common distance is supplied,
    distance = distance * ones(numel(az),1);        % use the same distance for every point.
else                                                % Otherwise,
    distance = distance(:);                         % treat distance as one value per point.
end

pos = zeros(numel(az),3);                           % Allocate output coordinates.
for ii = 1:numel(az)                                % Generate one point at a time.
    dir = [cos(el(ii))*cos(az(ii)), ...             % x component of unit direction.
           cos(el(ii))*sin(az(ii)), ...             % y component of unit direction.
           sin(el(ii))];                            % z component of unit direction.
    pos(ii,:) = center + distance(ii)*dir;          % Point coordinate with exact distance.
end
end
