# vim: set fileencoding=latin-1
#Copyright (c) 2003, 2008 Janne Blomqvist

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

from pylab import *


class Cell:
    """Class for representing a supercell."""
    
    def __init__(self, poscar=None, xyz=None):
        """Initialize the data members of this class"""
        # List of the chemical symbols of the atoms
        self.atomNames = []
        # Lattice constant, in Ångströms
        self.latticeConstant = 1.
        # 3x3 matrix containing the basis vectors of the supercell
        # in row major format
        self.basisVectors = eye(3)
        # Array containing the numbers of each element in the
        # system, i.e. the length of this array is the same as the
        # length of the self.atomNames list
        self.nAtomsType = array(0)
        # Are the ions allowed to move?
        self.selectiveDynamics = True
        # Flags for each atom describing in which cartesian coordinate
        # direction the atom is allowed to move. It is thus a 3xnAtoms
        # size list
        self.selectiveFlags = []
        # Are the atomic coordinates cartesian or in direct coordinates
        # If direct, cartesian coordinates can be calculated by
        # multiplying each coordinate with the basis vector matrix
        self.cartesian = True
        # Are the atomic coordinates relative or absolute, i.e. should
        # they be multiplies by the lattice constant.
        self.relative = False
        # Coordinates of the atoms
        self.atoms = zeros((0, 3))
        if (poscar != None):
            read_poscar(poscar)
        elif (xyz != None):
            read_xyz(xyz)

    def read_poscar(self, filename):
        """Parses a POSCAR file"""
        f = open(filename)
        poscar = f.readlines()
        f.close()
            
        # First line should contain the atom names , eg. "Ag Ge" in
        # the same order
        # as later in the file (and POTCAR for the full vasp run)
        self.atomNames = string.split(poscar[0])
            
        self.latticeConstant = float(poscar[1])
            
        # Now the lattice vectors
        a = []
        for vector in poscar[2:5]:
            s = string.split(vector)
            floatvect = float(s[0]), float(s[1]), float(s[2])
            a.append( floatvect)
        
        # Transpose to make natural ordering for linear algebra
        self.basisVectors = transpose(array(a))
        
        # Number of atoms. Again this must be in the same order as
        # in the first line
        # and in the POTCAR file
        numofatoms = string.split(poscar[5])
        tot_numatoms = 0
        for i in xrange(len(numofatoms)):
            numofatoms[i] = int(numofatoms[i])
        self.nAtomsType = array(numofatoms)
        self.nAtoms = sum(self.nAtomsType)
        
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
        atomcoords = []
        self.selectiveFlags = []
        for natomType in numofatoms:
            for atype in xrange(natomType):
                ac = string.split(poscar[atype+offset])
                atomcoords.append((float(ac[0]), float(ac[1]), float(ac[2])))
                if self.selectiveDynamics:
                    self.selectiveFlags.append((ac[3], ac[4], ac[5]))
            offset = offset + natomType
        
        # Transpose to produce sensible linear algebra
        self.atoms = transpose(array(atomcoords))

    def write_poscar(self, atoms, filename="POSCAR.out"):
        """Writes data into a POSCAR format file"""
        fc = "" # Contents of the file
        for a in atoms.atomNames:
            fc += str(a) + " "
        fc += "\n" + str(atoms.latticeConstant) + "\n"
        for i in range(0,3):
            for j in range(0,3):
                fc += str(atoms.basisVectors[i,j]) + " "
            fc += "\n"
        for at in atoms.nAtomsType:
            fc += str(at) + " "
        fc += "\n"
        if atoms.selectiveDynamics:
            fc += "Selective dynamics\n"
        if atoms.cartesian:
            fc += "Cartesian\n"
        else:
            fc += "Direct\n"
        for i in arange(0,atoms.nAtoms):
            for j in range(0,3):
                fc += str(atoms.atoms[j,i]) + " "
            if atoms.selectiveDynamics:
                selflags = atoms.selectiveFlags[i]
                for j in range(0,3):
                    fc += str(selflags[j]) + " "
            fc += "\n"
        f = open(filename, "w")
        f.write(fc)
        f.close()

    def write_xyz(self, atoms, filename="Xyz.out"):
        """Writes data into a XYZ format file"""
        fc = "" # Contents of the file
        fc += str(atoms.nAtoms) + "\nGenerated by vasputil\n"
        aindex = 0
        anameindex = 0
        if not atoms.cartesian:
            atoms.direct2Cartesian()
        for nAtomType in atoms.nAtomsType:
            for nAtom in arange(nAtomType):
                fc += atoms.atomNames[anameindex] + "\t"
                for i in range(3):
                    fc += str(atoms.atoms[i,aindex]) + "\t"
                fc += "\n"
                aindex += 1
            anameindex += 1

        f = open(filename, "w")
        f.write(fc)
        f.close()



    def cartesian2Direct(self):
        """Convert atom coordinates from cartesian to direct"""
        if not self.cartesian:
            return
        self.atoms = linalg.solve(self.latticeConstant*self.basisVectors, self.atoms)
        self.cartesian = False

    def direct2Cartesian(self):
        """Convert atom coordinates from direct to cartesian"""
        if self.cartesian:
            return
        self.atoms = dot(self.latticeConstant*self.basisVectors, self.atoms)
        self.cartesian = True
        
