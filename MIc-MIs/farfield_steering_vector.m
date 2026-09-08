function a = farfield_steering_vector(arrayPos, pointPos, lambda, normalizeFlag)
%FARFIELD_STEERING_VECTOR UPA/array steering vector from coordinates.
% arrayPos is an N-by-3 array of element coordinates. pointPos is a 1-by-3
% far-field point. The steering phase is generated from the propagation
% direction between the array center and the point. This implements the
% far-field UPA steering-vector idea in the paper while allowing arbitrary
% planar array coordinates.

if nargin < 4                                      % If normalization is not specified,
    normalizeFlag = false;                         % keep unit-modulus entries by default.
end
center = mean(arrayPos,1);                         % Array center.
linkVec = pointPos(:).' - center;                  % Direction vector from array center to point.
d0 = norm(linkVec);                                % Center-to-point distance.
u = linkVec / max(d0,1e-15);                       % Unit propagation direction.
relPos = arrayPos - center;                        % Element coordinates relative to array center.
a = exp(-1j*2*pi/lambda * (relPos*u.'));           % Far-field steering vector with unit-modulus entries.
if normalizeFlag                                   % If requested by the target-response model,
    a = a / max(norm(a),1e-15);                    % normalize the vector to unit norm.
end
end
