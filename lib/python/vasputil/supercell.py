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

import pylab as m


class Cell(object):
    """Class for representing a supercell."""
    
    def __init__(self, poscar=None, xyz=None):
        """Initialize the data members of this class"""
        # List of the chemical symbols of the atoms
        self.atomSymbols = []
        # Lattice constant, in Ångströms
        self.latticeConstant = 1.
        # 3x3 matrix containing the basis vectors of the supercell
        # in row major format
        self.basisVectors = m.eye(3)
        # Are the ions allowed to move?
        self.selectiveDynamics = False
        # Flags for each atom describing in which cartesian coordinate
        # direction the atom is allowed to move. It is thus a nAtomsx3
        # size list
        self.selectiveFlags = []
        # Are the atomic coordinates cartesian or in direct coordinates
        # If direct, cartesian coordinates can be calculated by
        # multiplying each coordinate with the basis vector matrix
        # times the lattice constant
        self.cartesian = True
        # Coordinates of the atoms
        self.atoms = m.zeros((0, 3))
        if (poscar != None):
            self.read_poscar(poscar)
        elif (xyz != None):
            self.read_xyz(xyz)

    def __getNAtoms(self):
        return self.atoms.shape[0]

    def __setNAtoms(self, val):
        raise AttributeError, "can't set attribute"

    def __delNAtoms(self):
        raise AttributeError, "can't delete attribute"

    nAtoms = property(__getNAtoms, __setNAtoms, __delNAtoms, "Total number of atoms.")

    def read_poscar(self, filename):
        """Parses a POSCAR file"""
        f = open(filename)
        poscar = f.readlines()
        f.close()
            
        # First line should contain the atom names , eg. "Ag Ge" in
        # the same order
        # as later in the file (and POTCAR for the full vasp run)
        atomNames = poscar[0].split()
            
        self.latticeConstant = float(poscar[1])
            
        # Now the lattice vectors
        a = []
        for vector in poscar[2:5]:
            s = vector.split()
            floatvect = float(s[0]), float(s[1]), float(s[2])
            a.append( floatvect)
        
        # Transpose to make natural ordering for linear algebra
        self.basisVectors = m.transpose(m.array(a))
        
        # Number of atoms. Again this must be in the same order as
        # in the first line
        # and in the POTCAR file
        numofatoms = poscar[5].split()
        for i in xrange(len(numofatoms)):
            numofatoms[i] = int(numofatoms[i])
            if (len(atomNames) < i + 1):
                atomNames.append("Unknown")
            [self.atomSymbols.append(atomNames[i]) for n in xrange(numofatoms[i])]
        
        # Check if Selective dynamics is switched on
        sdyn = poscar[6]
        add = 0
        if sdyn[0] == "S" or sdyn[0] == "s":
            add = 1
            self.selectiveDynamics = True
        
        # Check if atom coordinates are cartesian or direct
        acType = poscar[6+add]
        if acType[0] == "C" or acType[0] == "c" or acType[0] == "K" or acType[0] == "k":
            self.cartesian = 1
        else:
            self.cartesian = 0
        
        offset = add+7
        tot_natoms = sum(numofatoms)
        self.atoms = m.zeros((tot_natoms, 3))
        self.selectiveFlags = []
        for atom in xrange(tot_natoms):
            ac = poscar[atom+offset].split()
            self.atoms[atom] = (float(ac[0]), float(ac[1]), float(ac[2]))
            if self.selectiveDynamics:
                self.selectiveFlags.append((ac[3], ac[4], ac[5]))
        if self.cartesian:
            self.atoms *= self.latticeConstant
        

    def write_poscar(self, filename="POSCAR.out", fd=None):
        """Writes data into a POSCAR format file"""
        fc = "" # Contents of the file
        asc = self.getAtomSymbolCount()
        for a in asc:
            fc += str(a[0]) + " "
        fc += "\n" + "%19.16f" % self.latticeConstant + "\n"
        basisfmt = "%22.16f"
        for i in xrange(3):
            fc += " "
            for j in xrange(3):
                fc += basisfmt % self.basisVectors[i,j] 
            fc += "\n"
        atomnumfmt = "%4i"
        for at in asc:
            fc += atomnumfmt % at[1]
        fc += "\n"
        if self.selectiveDynamics:
            fc += "Selective dynamics\n"
        if self.cartesian:
            fc += "Cartesian\n"
            atoms = self.atoms / self.latticeConstant
        else:
            fc += "Direct\n"
            atoms = self.atoms
        atomfmt = "%20.16f"
        logicalfmt = "%4s"
        for i in xrange(self.nAtoms):
            for j in xrange(3):
                fc += atomfmt % atoms[i,j]
            if self.selectiveDynamics:
                selflags = self.selectiveFlags[i]
                for j in xrange(3):
                    fc += logicalfmt % selflags[j]
            fc += "\n"
        if (fd == None):
            f = open(filename, "w")
            f.write(fc)
            f.close()
        else:
            fd.write(fc)
        
    def read_xyz(self, infile):
        "Parses an xyz file"
        f = open(infile)
        xyz = f.readlines()
        f.close()
        # first line contains number of atoms
        self.atoms = m.zeros((int(xyz[0]), 3))
        self.cartesian = True
        self.atomSymbols = []
        for ii in xrange(2, self.nAtoms):
            s = xyz[ii].split()
            floatvect = m.array([float(s[1]), float(s[2]), float(s[3])])
            self.atoms[(ii-1),:] = floatvect
            self.atomSymbols.append(s[0])
        return self.atoms


    def write_xyz(self, filename="atoms.xyz", comment="Generated by vasputil"):
        """Writes data into a XYZ format file.
        
        Format is the same as used by OpenBabel, see
        http://openbabel.org/wiki/XYZ_(format)
        
        """
        fc = "" # Contents of the file
        fc += str(self.nAtoms) + "\n" + comment + "\n"
        if not self.cartesian:
            self.direct2Cartesian()
        for nn in xrange(self.nAtoms):
            fc += "%3s" % self.atomSymbols[nn]
            for i in xrange(3):
                fc += "%10.5f" % self.atoms[nn,i]
            fc += "\n"
        f = open(filename, "w")
        f.write(fc)
        f.close()

    def sortAtoms(self):
        """Sort the atoms in the cell according to atomic symbol."""
        ind = m.argsort(self.atomSymbols)
        self.atomSymbols = m.array(self.atomSymbols)[ind]
        self.atoms = self.atoms[ind]
        self.selectiveFlags = m.array(self.selectiveFlags)[ind]

    def getAtomSymbolCount(self):
        """Return a list of (atomic symbol, count) tuples.

        As a side effect, sorts the atoms by symbol.

        """
        self.sortAtoms()
        sc = []
        psym = self.atomSymbols[0]
        count = 0
        for sym in self.atomSymbols:
            if (sym != psym):
                sc.append((psym, count))
                psym = sym
                count = 1
            else:
                count += 1
        sc.append((psym, count))
        return sc

    def cartesian2Direct(self):
        """Convert atom coordinates from cartesian to direct"""
        if not self.cartesian:
            return
        self.atoms = m.transpose(m.linalg.solve(self.latticeConstant * \
                self.basisVectors, \
                m.transpose(self.atoms)))
        self.cartesian = False

    def direct2Cartesian(self):
        """Convert atom coordinates from direct to cartesian"""
        if self.cartesian:
            return
        self.atoms = m.transpose(m.dot(self.latticeConstant*self.basisVectors, \
                m.transpose(self.atoms)))
        self.cartesian = True
        
    def showVmd(self):
        """Show a supercell in VMD."""
# This is a quick and dirty hack, as VMD has some builtin support 
# as well.
        import tempfile, os
        f = tempfile.NamedTemporaryFile()
        self.write_poscar(fd=f)
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
        self.direct2Cartesian()
        cell.direct2Cartesian()
        oldsz = self.nAtoms
        atoms = m.zeros((self.nAtoms + cell.nAtoms, 3))
        self.atomSymbols = m.array(self.atomSymbols + cell.atomSymbols)
        self.selectiveFlags = m.array(self.selectiveFlags + cell.selectiveFlags)
        atoms[:oldsz] = self.atoms
        atoms[oldsz:] = cell.atoms
        self.atoms = atoms


# End of class Cell

def rotate_molecule(coords, rotp = m.array((0.,0.,0.)), phi = 0., \
        theta = 0., psi = 0.):
    """Rotate a molecule via Euler angles.

    See http://mathworld.wolfram.com/EulerAngles.html for definition.
    Input arguments:
    coords: Atom coordinates, as Nx3 2d pylab array.
    rotp: The point to rotate about, as a 1d 3-element pylab array
    phi: The 1st rotation angle around z axis.
    theta: Rotation around x axis.
    psi: 2nd rotation around z axis.

    """
# First move the molecule to the origin
# In contrast to MATLAB, numpy broadcasts the smaller array to the larger
# row-wise, so there is no need to play with the Kronecker product.
    rcoords = coords - rotp
# First Euler rotation about z in matrix form
    D = m.array(((m.cos(phi), m.sin(phi), 0.), (-m.sin(phi), m.cos(phi), 0.), \
            (0., 0., 1.)))
# Second Euler rotation about x:
    C = m.array(((1., 0., 0.), (0., m.cos(theta), m.sin(theta)), \
            (0., -m.sin(theta), m.cos(theta))))
# Third Euler rotation, 2nd rotation about z:
    B = m.array(((m.cos(psi), m.sin(psi), 0.), (-m.sin(psi), m.cos(psi), 0.), \
            (0., 0., 1.)))
# Total Euler rotation
    A = m.dot(B, m.dot(C, D))
# Do the rotation
    rcoords = m.dot(A, m.transpose(rcoords))
# Move back to the rotation point
    return m.transpose(rcoords) + rotp
