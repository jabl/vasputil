# vim: set fileencoding=latin-1
# Copyright (c) 2003, 2008 Janne Blomqvist

#  This file is part of Vasputil.

#  Vasputil is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.

#  Vasputil is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with vasputil.  If not, see <http://www.gnu.org/licenses/>.

"""This module defines a utility functions for dealing with vasp 
supercells. Instead of the old Cell class, the module now uses
the Atoms class from ase. 

"""

try:
    import numpy as np
except ImportError:
    import pylab as np

import vasputil.geometry as vg

def natoms(atoms):
    """How many atoms in the cell?"""
    return atoms.positions.shape[0]

def atoms_distance(atoms, atom1, atom2, proj=None):
    """Measure the distance between two atoms.
    
    Atoms are indexed starting from 0, following the usual Python
    convention.  Note that this is different from VASP itself, which starts
    indexing from 1.  This method takes into account periodic boundary
    conditions.

    Arguments:
    atoms -- ASE Atoms object containing the atoms.
    atom1 -- The index of one of the atoms, starting from 0.
    atom2 -- The index of the other atom, starting from 0.
    proj  -- Projection along a vector or plane. If a string, it can
             contain x, y, z and the method then measures the distance
             in the plane defined by the string. If it's a sequence
             of three numbers, the method measures the distance
             projected along the vector.
    
    """
    at = atoms.get_scaled_positions()
    dvec = at[atom1, :] - at[atom2, :]
    dvec = np.dot(vg.vec_pbc(dvec), \
            atoms.get_cell())
    if proj == None:
        return np.linalg.norm(dvec)
    elif type(proj) == str:
        if len(proj) != 2:
            raise TypeError("Length of string specifying plane must be 2.")
        pvec = dvec.copy()
        if proj.find("x") == -1:
            pvec[0] = 0.
        if proj.find("y") == -1:
            pvec[1] = 0.
        if proj.find("z") == -1:
            pvec[2] = 0.
        return abs(np.dot(dvec, pvec) / np.linalg.norm(pvec))
    else:
        return abs(np.dot(dvec, proj) / np.linalg.norm(proj))

def nearest_neighbors(atoms, tol=1.0, num_neigh=None):
    """Nearest neighbors and distances.

    Arguments:
    atoms -- The ASE Atoms object with all the atoms.
    tol  -- Return only distances smaller than this. Default 1.0 Å.
    num_neigh -- Number of nearest neighbors per atom returned.

    Returns -- List containing 
               (source_atom, target_atom, dist) tuples. 

    """
    at = atoms.get_scaled_positions()
    nn = []
    for anum in range(len(at)):
        dvec = at - at[anum]
        dvec = np.dot(vg.vec_pbc(dvec), \
                atoms.get_cell())
        dist = np.empty(dvec.shape[0])
        for ii in range(len(dvec)):
            dist[ii] = np.linalg.norm(dvec[ii])
        if num_neigh == None:
            mask = dist < tol
            for ii in range(len(mask)):
                if mask[ii] and ii != anum:
                    nn.append((anum, ii, dist[ii]))
        else:
            sind = dist.argsort()
            for ii in range(min(num_neigh + 1, len(dist))):
                if anum != sind[ii]:
                    nn.append((anum, sind[ii], dist[sind[ii]]))
    return nn

def atoms_moved(cell1, cell2, tol=0.1):
    """Return a list of atoms that have moved between the two cells.

    If lattices are compatible, take periodic boundary conditions into account.
    
    Arguments:
    cell1,2 -- The supercells to compare
    tol -- The tolerance in Å

    Return value -- A list of (atom index, distance moved) tuples.

    """
    (latt, nat) = check_cells(cell1, cell2)
    if latt:
        at1 = cell1.get_scaled_positions()
        at2 = cell2.get_scaled_positions()
    else:
        at1 = cell1.positions
        at2 = cell2.positions
    nmax = min(natoms(cell1), natoms(cell2))
    am = []
    for nn in range(nmax):
        dvec = at1[nn, :] - at2[nn, :]
        if latt:
            dvec = np.dot(cell1.get_cell(), \
                    vg.vec_pbc(dvec))
        dist = np.linalg.norm(dvec)
        if dist > tol:
            am.append((nn, dist))
    return am

def check_cells(cell1, cell2):
    """Check to which extent two cells are compatible.
    
    Return value -- a tuple where the first element is a boolean specifying
    whether the lattices are compatible, that is, comparing the basis vectors *
    lattice constants. The second element is a boolean specifying whether the
    cells contain an equal amount of atoms.
    
    """
    # First check that lattice constant * basis vectors are compatible.
    latt = np.any(cell1.get_cell() \
            - cell2.get_cell() < 1e-15)
    # Then check that there are an equal number of atoms.
    nat = natoms(cell1) == natoms(cell2)
    return (latt, nat)

def interpolate_cells(cell1, cell2, frac=0.5, images=1):
    """Interpolate coordinates between two supercells.
    
    Arguments:
    cell1 -- The starting point cell.
    cell2 -- The endpoint cell.
    frac -- Fraction, where on the interval [cell1,cell2] should the new cell
            reside. If 0.0, the resulting cell is equal to cell1, if 1.0 it's 
            equal to cell2.
    images -- Number of intermediate images. If != 1, frac is ignored.

    Return value -- A new cell with the interpolated coordinates, or a list 
                    of cells if images != 1.
    
    """
    import copy
    (latt, atoms) = check_cells(cell1, cell2)
    if not latt or not atoms:
        raise Error("Cells are not compatible.")
    if images == 1:
        icell = copy.deepcopy(cell1)
        icell.set_scaled_positions((1 - frac) * cell1.get_scaled_positions() \
                + frac * cell2.get_scaled_positions())
        return icell
    icells = []
    images += 1
    for ii in range(1, images):
        icell = copy.deepcopy(cell1)
        fr = float(ii) / images 
        icell.set_scaled_positions((1 - fr) * cell1.get_scaled_positions() \
                + fr * cell2.get_scaled_positions())
        icells.append(icell)
    return icells

def rotate_molecule(coords, rotp = np.array((0.,0.,0.)), phi = 0., \
        theta = 0., psi = 0.):
    """Rotate a molecule via Euler angles.

    See http://mathworld.wolfram.com/EulerAngles.html for definition.
    Input arguments:
    coords: Atom coordinates, as Nx3 2d numpy array.
    rotp: The point to rotate about, as a 1d 3-element numpy array
    phi: The 1st rotation angle around z axis.
    theta: Rotation around x axis.
    psi: 2nd rotation around z axis.

    """
    # First move the molecule to the origin
    # In contrast to MATLAB, numpy broadcasts the smaller array to the larger
    # row-wise, so there is no need to play with the Kronecker product.
    rcoords = coords - rotp
    # First Euler rotation about z in matrix form
    D = np.array(((np.cos(phi), np.sin(phi), 0.), (-np.sin(phi), np.cos(phi), 0.), \
            (0., 0., 1.)))
    # Second Euler rotation about x:
    C = np.array(((1., 0., 0.), (0., np.cos(theta), np.sin(theta)), \
            (0., -np.sin(theta), np.cos(theta))))
    # Third Euler rotation, 2nd rotation about z:
    B = np.array(((np.cos(psi), np.sin(psi), 0.), (-np.sin(psi), np.cos(psi), 0.), \
            (0., 0., 1.)))
    # Total Euler rotation
    A = np.dot(B, np.dot(C, D))
    # Do the rotation
    rcoords = np.dot(A, np.transpose(rcoords))
    # Move back to the rotation point
    return np.transpose(rcoords) + rotp
