function pos = upa_positions(M, spacing, center, plane)
%UPA_POSITIONS Generate square UPA coordinates.

side = round(sqrt(M));
if side^2 ~= M
    error('M must be a square number in this demo.');
end
grid = ((0:side-1) - (side-1)/2)*spacing;
pos = zeros(M,3);
idx = 0;
for iy = 1:side
    for iz = 1:side
        idx = idx + 1;
        switch lower(plane)
            case 'yz'
                pos(idx,:) = center + [0, grid(iy), grid(iz)];
            otherwise
                error('Only yz plane is implemented.');
        end
    end
end
end
