function rRayleigh = rayleigh_distance(pos, lambda)
%RAYLEIGH_DISTANCE Far-field boundary 2D^2/lambda for an array aperture.

D = aperture_diameter(pos);                         % Largest aperture dimension.
rRayleigh = 2*D^2/lambda;                           % Rayleigh far-field distance.
end
