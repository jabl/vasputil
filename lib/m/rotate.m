% rotate (coords, rotp, phi, theta, psi)
% Rotate molecule 
% See http://mathworld.wolfram.com/EulerAngles.html for definition
% of angles.
% Inputs:
% coords: The coordinates of the atoms in the molecule, as column
% vectors in a matrix.
% rotp: The point to rotate about. 
% phi: The 1st rotation angle around z axis.
% theta: Rotation around x axis.
% psi: 2nd rotation around z axis.
function rcoords = rotate (coords, rotp, phi, theta, psi)

% rotation point appropriately dimensioned:
dr = kron (ones (1, size (coords', 2)), rotp');

% First move the molecule by -rotp.
rcoords = coords' - dr;

% The first Euler rotation about z in matrix form:
D = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];

% The second Euler rotation, about x:
C = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];

% The third Euler rotation, about z again:
B = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];

A = B*C*D;

% Do the rotation
rcoords = A*rcoords;

% Move back to the rotation point
rcoords = rcoords + dr;
rcoords = rcoords';
