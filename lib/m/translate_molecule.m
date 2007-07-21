% translate_molecule (coords, dr)
% Translate a set of coordinates by a vector
% Inputs:
% coords: The coordinates, as row vectors.
% r: The translation vector
function tr = translate_molecule (coords, r)

% Rotation vector appropriately dimensioned
dr = kron (ones (size (coords, 1),1), r);

% Translate
tr = coords + dr;
