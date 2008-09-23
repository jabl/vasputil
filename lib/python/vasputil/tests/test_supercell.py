# vim: set fileencoding=latin-1
# Copyright (c) 2008 Janne Blomqvist

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

"""This module contains unit tests for the vasputil.supercell module."""

import unittest
import vasputil.supercell as s
import numpy.testing.numpytest as nt
import numpy.testing.utils as ntu
import math
import os

def testdir():
    """The directory where the tests and data files reside."""
    # WARNING: Ugly ugly ugly!
    # find the path where this file is
    return os.path.split(__file__)[0]

class CellTestCase(unittest.TestCase):
    """Testcase for vasputil.supercell.Cell class."""

    def setUp(self):
        # The test POSCAR file is in the same directory.
        path = os.path.join(testdir(), "POSCAR")
        self.cell = s.Cell(poscar=path)

    def test_lc(self):
        """Check that the lattice constant is imported correctly."""
        self.assertEqual(self.cell.lattice_constant, 1.)

    def test_atoms_distance(self):
        """Test the atoms_distance method."""
        dist = self.cell.atoms_distance(9, 24)
        ntu.assert_almost_equal(dist, 1.90823040889809)

    def test_atoms_distance_proj(self):
        """Test the atoms_distance method with a projection."""
        dist = self.cell.atoms_distance(24, 9, "xy")
        ntu.assert_almost_equal(dist, 1.51230705052) 

    def test_atoms_distance_proj2(self):
        """Test the atoms_distance method with a projection defined by vector."""
        dist = self.cell.atoms_distance(9, 24, (0,0,10))
        ntu.assert_almost_equal(dist, 1.16373136) 

    def test_atoms_distance_pbc(self):
        dist = self.cell.atoms_distance(2, 18)
        ntu.assert_almost_equal(dist, 1.867066, decimal=5)

class RotateMolTestCase(nt.NumpyTestCase):
    """Test the rotate_molecule function."""

    def test_rotate(self):
        coords = ntu.rand(10, 3)
        rcoords = s.rotate_molecule(coords, phi=math.pi, theta=2*math.pi, \
                psi=math.pi)
        ntu.assert_array_almost_equal(coords, rcoords)

class AtomsMovedTestCase(unittest.TestCase):
    """Test the atoms_moved function."""

    def setUp(self):
        p1 = os.path.join(testdir(), "POSCAR")
        self.c1 = s.Cell(poscar=p1)
        p2 = os.path.join(testdir(), "POSCAR2")
        self.c2 = s.Cell(poscar=p2)

    def test_moved(self):
        atoms = []
        for atom in s.atoms_moved(self.c1, self.c2):
            atoms.append(atom[0])
        self.failUnless(29 in atoms)

    def test_moved_tol(self):
        atoms = s.atoms_moved(self.c1, self.c2, 3)
        self.failUnless(len(atoms) == 0)

class CheckCellsTestCase(unittest.TestCase):
    """Test the check_cells function."""

    def setUp(self):
        p1 = os.path.join(testdir(), "POSCAR")
        self.c1 = s.Cell(poscar=p1)
        p2 = os.path.join(testdir(), "POSCAR2")
        self.c2 = s.Cell(poscar=p2)

    def test_check_cells(self):
        (latt, natoms) = s.check_cells(self.c1, self.c2)
        self.assertEqual(latt, True)
        self.assertEqual(natoms, True)

def suite():
    cell_suite = unittest.TestLoader().loadTestsFromTestCase(CellTestCase)
    rot_suite = unittest.TestLoader().loadTestsFromTestCase(RotateMolTestCase)
    move_suite = unittest.TestLoader().loadTestsFromTestCase(\
            AtomsMovedTestCase)
    check_suite = unittest.TestLoader().loadTestsFromTestCase(\
            CheckCellsTestCase)
    return unittest.TestSuite([cell_suite, rot_suite, move_suite, \
            check_suite])


if __name__ == "__main__":
    unittest.main()
