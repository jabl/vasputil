% Rotate molecule about z axis (it should be easy to generalize
% this to the full euler angle rotations).
% Inputs:
% coords: The coordinates of the atoms in the molecule, as column
% vectors in a matrix.
% rotp: The point to rotate about. Currently only the x and y
% coordinate are used.
% angle: the rotation angle.
function rcoords = rotatez (coords, rotp, angle)

% rotation point appropriately dimensioned:
dr = kron (ones (1, size (coords, 2)), rotp);

% First move the molecule by -rotp.
rcoords = coords - dr;

% The first Euler rotation in matrix form:
A = [cos(angle) sin(angle) 0; -sin(angle) cos(angle) 0; 0 0 1];

% Do the rotation
rcoords = A*rcoords;

% Move back to the rotation point
rcoords = rcoords + dr;
