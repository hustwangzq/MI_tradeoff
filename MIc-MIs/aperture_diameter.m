function D = aperture_diameter(pos)
%APERTURE_DIAMETER Maximum linear aperture size of a planar array.
% For an Nr-by-Nc UPA with half-wavelength spacing, this gives the largest
% edge length rather than the diagonal.  The Rayleigh distance is 2D^2/lambda.

spanX = max(pos(:,1))-min(pos(:,1));                % Array span along x.
spanY = max(pos(:,2))-min(pos(:,2));                % Array span along y.
spanZ = max(pos(:,3))-min(pos(:,3));                % Array span along z.
D = max([spanX, spanY, spanZ]);                     % Largest physical aperture dimension.
end
