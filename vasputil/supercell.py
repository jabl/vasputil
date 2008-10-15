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

"""This module defines a class that represents a supercell, as well as utility
functions.

"""

try:
    import numpy as n
except ImportError:
    import pylab as n

import vasputil.geometry as vg


class Cell(object):
    """Class for representing a supercell."""
    
    def __init__(self, poscar=None, xyz=None):
        """Initialize the data members of this class"""
        # List of the chemical symbols of the atoms
        self.atom_symbols = []
        # Lattice constant, in Ångströms
        self.lattice_constant = 1.
        # 3x3 matrix containing the basis vectors of the supercell
        # as row vectors
        self.basis_vectors = n.eye(3)
        # Are the ions allowed to move?
        self.selective_dynamics = False
        # Flags for each atom describing in which cartesian coordinate
        # direction the atom is allowed to move. It is thus a natomsx3
        # size list
        self.selective_flags = []
        # Are the atomic coordinates cartesian or in direct coordinates
        # If direct, cartesian coordinates can be calculated by
        # multiplying each coordinate with the basis vector matrix
        # times the lattice constant
        self.cartesian = True
        # Coordinates of the atoms
        self.atoms = n.zeros((0, 3))
        if (poscar != None):
            self.read_poscar(poscar)
        elif (xyz != None):
            self.read_xyz(xyz)

    def __get_natoms(self):
        return self.atoms.shape[0]

    def __set_natoms(self, val):
        raise AttributeError, "can't set attribute"

    def __del_natoms(self):
        raise AttributeError, "can't delete attribute"

    natoms = property(__get_natoms, __set_natoms, __del_natoms, \
            "Total number of atoms.")

    def _get_volume(self):
        return n.dot(self.basis_vectors[:,0], n.cross( \
                self.basis_vectors[:,1], self.basis_vectors[:,2]))

    def _set_volume(self, val):
        raise AttributeError("Can't set attribute.")

    def _del_volume(self):
        raise AttributeError("Can't delete attribute.")

    volume = property(_get_volume, _set_volume, _del_volume, \
            "The volume of the supercell.")

    def read_poscar(self, pfile):
        """Parses a POSCAR file"""
        if type(pfile) == str:
            f = open(pfile)
        elif type(pfile) == file:
            f = pfile
        else:
            raise TypeError("pfile argument must be a string or a file object.")
            
        # First line should contain the atom names , eg. "Ag Ge" in
        # the same order
        # as later in the file (and POTCAR for the full vasp run)
        atomNames = f.readline().split()
            
        self.lattice_constant = float(f.readline())
            
        # Now the lattice vectors
        a = []
        for ii in range(3):
            s = f.readline().split()
            floatvect = float(s[0]), float(s[1]), float(s[2])
            a.append(floatvect)
        
        self.basis_vectors = n.array(a)
        
        # Number of atoms. Again this must be in the same order as
        # in the first line
        # and in the POTCAR file
        numofatoms = f.readline().split()
        for i in xrange(len(numofatoms)):
            numofatoms[i] = int(numofatoms[i])
            if (len(atomNames) < i + 1):
                atomNames.append("Unknown")
            [self.atom_symbols.append(atomNames[i]) \
                    for na in xrange(numofatoms[i])]
        
        # Check if Selective dynamics is switched on
        sdyn = f.readline()
        if sdyn[0] == "S" or sdyn[0] == "s":
            self.selective_dynamics = True
        
        # Check if atom coordinates are cartesian or direct
        if self.selective_dynamics:
            ac_type = f.readline()
        else:
            ac_type = sdyn
        if ac_type[0] == "C" or ac_type[0] == "c" or ac_type[0] == "K" \
                or ac_type[0] == "k":
            self.cartesian = 1
        else:
            self.cartesian = 0
        tot_natoms = sum(numofatoms)
        self.atoms = n.empty((tot_natoms, 3))
        self.selective_flags = []
        for atom in xrange(tot_natoms):
            ac = f.readline().split()
            self.atoms[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
            if self.selective_dynamics:
                self.selective_flags.append((ac[3], ac[4], ac[5]))
        if self.cartesian:
            self.atoms *= self.lattice_constant
        if type(pfile) == str:
            f.close()
        

    def write_poscar(self, pfile):
        """Writes data into a POSCAR format file"""
        fc = "" # Contents of the file
        asc = self.get_atom_symbol_count()
        for a in asc:
            fc += str(a[0]) + " "
        fc += "\n" + "%19.16f" % self.lattice_constant + "\n"
        basisfmt = "%22.16f"
        for i in xrange(3):
            fc += " "
            for j in xrange(3):
                fc += basisfmt % self.basis_vectors[i,j] 
            fc += "\n"
        atomnumfmt = "%4i"
        for at in asc:
            fc += atomnumfmt % at[1]
        fc += "\n"
        if self.selective_dynamics:
            fc += "Selective dynamics\n"
        if self.cartesian:
            fc += "Cartesian\n"
            atoms = self.atoms / self.lattice_constant
        else:
            fc += "Direct\n"
            atoms = self.atoms
        atomfmt = "%20.16f"
        logicalfmt = "%4s"
        for i in xrange(self.natoms):
            for j in xrange(3):
                fc += atomfmt % atoms[i,j]
            if self.selective_dynamics:
                selflags = self.selective_flags[i]
                for j in xrange(3):
                    fc += logicalfmt % selflags[j]
            fc += "\n"
        if type(pfile) == str:
            f = open(pfile, "w")
            f.write(fc)
            f.close()
        elif type(pfile) == file:
            pfile.write(fc)
        else:
            raise TypeError("pfile argument must be string of file object.")
        
    def read_xyz(self, infile):
        "Parses an xyz file"
        f = open(infile)
        xyz = f.readlines()
        f.close()
        # first line contains number of atoms
        self.atoms = n.empty((int(xyz[0]), 3))
        self.cartesian = True
        self.atom_symbols = []
        for ii in xrange(2, self.natoms):
            s = xyz[ii].split()
            floatvect = n.array([float(s[1]), float(s[2]), float(s[3])])
            self.atoms[(ii-1),:] = floatvect
            self.atom_symbols.append(s[0])
        return self.atoms


    def write_xyz(self, filename="atoms.xyz", comment="Generated by vasputil"):
        """Writes data into a XYZ format file.
        
        Format is the same as used by OpenBabel, see
        http://openbabel.org/wiki/XYZ_(format)
        
        """
        fc = "" # Contents of the file
        fc += str(self.natoms) + "\n" + comment + "\n"
        self.direct2cartesian()
        for nn in xrange(self.natoms):
            fc += "%3s" % self.atom_symbols[nn]
            for i in xrange(3):
                fc += "%10.5f" % self.atoms[nn,i]
            fc += "\n"
        f = open(filename, "w")
        f.write(fc)
        f.close()

    def sort_atoms(self):
        """Sort the atoms in the cell according to atomic symbol."""
        ind = n.argsort(self.atom_symbols)
        self.atom_symbols = n.array(self.atom_symbols)[ind]
        self.atoms = self.atoms[ind]
        self.selective_flags = n.array(self.selective_flags)[ind]

    def get_atom_symbol_count(self):
        """Return a list of (atomic symbol, count) tuples.

        As a side effect, sorts the atoms by symbol.

        """
        self.sort_atoms()
        sc = []
        psym = self.atom_symbols[0]
        count = 0
        for sym in self.atom_symbols:
            if (sym != psym):
                sc.append((psym, count))
                psym = sym
                count = 1
            else:
                count += 1
        sc.append((psym, count))
        return sc

    def cartesian2direct(self):
        """Convert atom coordinates from cartesian to direct.
        
        For further explanation, see documentation for direct2cartesian.
        In this case, we want to solve atoms_d * basis = atoms_c for atoms_d.
        In order to turn it into the canonical Ax = b system we need to
        transpose everything, then solve, and transpose back.
        
        """
        if not self.cartesian:
            return
        self.atoms = n.transpose(n.linalg.solve(self.lattice_constant * \
                n.transpose(self.basis_vectors), \
                n.transpose(self.atoms)))
        self.cartesian = False

    def direct2cartesian(self):
        """Convert atom coordinates from direct to cartesian.
        
        If the basis vectors are stored as column vectors in an array, the
        coordinates of an atom as a column vector can be converted from
        direct to cartesian coordinates by basis * atom_coords. Generalizing
        to multiple atoms implies that atom_coords must then be a 3xN array.

        Now, we store atoms as a Nx3 array, and basis vectors as row vectors,
        so both of these need to be transposed. However, from linear algebra
        we know that (basis.T * atoms.T).T = atoms * basis, where .T is the
        transpose operator.
        
        """
        if self.cartesian:
            return
        self.atoms = n.dot(self.atoms, \
                self.lattice_constant * self.basis_vectors)
        self.cartesian = True
        
    def show_vmd(self):
        """Show a supercell in VMD."""
        # This is a quick and dirty hack, as VMD has some builtin support 
        # as well.
        import tempfile, os
        f = tempfile.NamedTemporaryFile()
        self.write_poscar(f)
        f.flush()
        vmdstr = "vmd -nt -POSCAR " + f.name
        print "now executing " + vmdstr
        os.system(vmdstr)
        raw_input("Press Enter when done to delete the temp file.")
        f.close()

    def add(self, cell):
        """Add another supercell to this one. 
        
        The resulting cell will have
        the basis vector and lattice constant of this cell.
        
        """
        self.direct2cartesian()
        cell.direct2cartesian()
        oldsz = self.natoms
        atoms = n.empty((self.natoms + cell.natoms, 3))
        self.atom_symbols = n.array(self.atom_symbols + cell.atom_symbols)
        self.selective_flags = n.array(self.selective_flags + cell.selective_flags)
        atoms[:oldsz] = self.atoms
        atoms[oldsz:] = cell.atoms
        self.atoms = atoms

    def atoms_distance(self, atom1, atom2, proj=None):
        """Measure the distance between two atoms.
        
        Atoms are indexed starting from 0, following the usual Python
        convention.  Note that this is different from VASP itself, which starts
        indexing from 1.  This method takes into account periodic boundary
        conditions.

        Arguments:
        atom1 -- The index of one of the atoms, starting from 0.
        atom2 -- The index of the other atom, starting from 0.
        proj  -- Projection along a vector or plane. If a string, it can
                 contain x, y, z and the method then measures the distance
                 in the plane defined by the string. If it's a sequence
                 of three numbers, the method measures the distance
                 projected along the vector.
        
        """
        self.cartesian2direct()
        dvec = self.atoms[atom1, :] - self.atoms[atom2, :]
        dvec = n.dot(vg.vec_pbc(dvec), \
                self.lattice_constant * self.basis_vectors)
        if proj == None:
            return n.linalg.norm(dvec)
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
            return abs(n.dot(dvec, pvec) / n.linalg.norm(pvec))
        else:
            return abs(n.dot(dvec, proj) / n.linalg.norm(proj))

    def nearest_neighbors(self, tol=1.0, num_neigh=None):
        """Nearest neighbors and distances.

        Arguments:
        tol  -- Return only distances smaller than this. Default 1.0 Å.
        num_neigh -- Number of nearest neighbors per atom returned.

        Returns -- List containing 
                   (source_atom, target_atom, dist) tuples. 

        """
        self.cartesian2direct()
        nn = []
        for anum in range(len(self.atoms)):
            dvec = self.atoms - self.atoms[anum]
            dvec = n.dot(vg.vec_pbc(dvec), \
                    self.lattice_constant * self.basis_vectors)
            dist = n.empty(dvec.shape[0])
            for ii in range(len(dvec)):
                dist[ii] = n.linalg.norm(dvec[ii])
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
        


# End of class Cell

def atoms_moved(cell1, cell2, tol=0.1):
    """Return a list of atoms that have moved between the two cells.

    If lattices are compatible, take periodic boundary conditions into account.
    
    Arguments:
    cell1,2 -- The supercells to compare
    tol -- The tolerance in Å

    Return value -- A list of (atom index, distance moved) tuples.

    """
    (latt, natoms) = check_cells(cell1, cell2)
    if latt:
        cell1.cartesian2direct()
        cell2.cartesian2direct()
    else:
        cell1.direct2cartesian()
        cell2.direct2cartesian()
    nmax = min(cell1.natoms, cell2.natoms)
    am = []
    for nn in range(nmax):
        dvec = cell1.atoms[nn, :] - cell2.atoms[nn, :]
        if latt:
            dvec = n.dot(cell1.lattice_constant * cell1.basis_vectors, \
                    vg.vec_pbc(dvec))
        dist = n.linalg.norm(dvec)
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
    latt = n.any(cell1.lattice_constant * cell1.basis_vectors \
            - cell2.lattice_constant * cell2.basis_vectors < 1e-15)
    # Then check that there are an equal number of atoms.
    natoms = cell1.natoms == cell2.natoms
    return (latt, natoms)

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
    cell1.cartesian2direct()
    cell2.cartesian2direct()
    if images == 1:
        icell = copy.copy(cell1)
        icell.atoms = (1 - frac) * cell1.atoms + frac * cell2.atoms
        return icell
    icells = []
    images += 1
    for ii in range(1, images):
        icell = copy.copy(cell1)
        fr = float(ii) / images 
        icell.atoms = (1 - fr) * cell1.atoms + fr * cell2.atoms
        icells.append(icell)
    return icells

def rotate_molecule(coords, rotp = n.array((0.,0.,0.)), phi = 0., \
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
    D = n.array(((n.cos(phi), n.sin(phi), 0.), (-n.sin(phi), n.cos(phi), 0.), \
            (0., 0., 1.)))
    # Second Euler rotation about x:
    C = n.array(((1., 0., 0.), (0., n.cos(theta), n.sin(theta)), \
            (0., -n.sin(theta), n.cos(theta))))
    # Third Euler rotation, 2nd rotation about z:
    B = n.array(((n.cos(psi), n.sin(psi), 0.), (-n.sin(psi), n.cos(psi), 0.), \
            (0., 0., 1.)))
    # Total Euler rotation
    A = n.dot(B, n.dot(C, D))
    # Do the rotation
    rcoords = n.dot(A, n.transpose(rcoords))
    # Move back to the rotation point
    return n.transpose(rcoords) + rotp
