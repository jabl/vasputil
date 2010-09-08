# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8

# Copyright (c) 2008, 2010 Janne Blomqvist

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
from ase.io import read as aread
import numpy as np
import math
import os

def testdir():
    """The directory where the tests and data files reside."""
    # WARNING: Ugly ugly ugly!
    # find the path where this file is
    return os.path.split(__file__)[0]

def load_cells():
    """Load the two test supercells."""
    p1 = os.path.join(testdir(), "POSCAR")
    c1 = aread(p1, format="vasp")
    p2 = os.path.join(testdir(), "POSCAR2")
    c2 = aread(p2, format="vasp")
    return (c1, c2)

def array_almost_equal(a1, a2, tol=np.finfo(type(1.0)).eps):
    """Replacement for old numpy.testing.utils.array_almost_equal."""
    return (np.abs(a1 - a2) < tol).all()

class CellTestCase(unittest.TestCase):
    """Testcase for vasputil.supercell.Cell class."""

    def setUp(self):
        # The test POSCAR file is in the same directory.
        path = os.path.join(testdir(), "POSCAR")
        self.cell = aread(path, format="vasp")

    def test_atoms_distance(self):
        """Test the atoms_distance method."""
        dist = s.atoms_distance(self.cell, 9, 24)
        self.assertAlmostEqual(dist, 1.90823040889809)

    def test_atoms_distance_proj(self):
        """Test the atoms_distance method with a projection."""
        dist = s.atoms_distance(self.cell, 24, 9, "xy")
        self.assertAlmostEqual(dist, 1.51230705052) 

    def test_atoms_distance_proj2(self):
        """Test the atoms_distance method with a projection defined by vector."""
        dist = s.atoms_distance(self.cell, 9, 24, (0,0,10))
        self.assertAlmostEqual(dist, 1.16373136) 

    def test_atoms_distance_pbc(self):
        dist = s.atoms_distance(self.cell, 2, 18)
        self.assertAlmostEqual(dist, 1.867066, 5)

    def test_nndist(self):
        nnd = s.nearest_neighbors(self.cell, tol=1.7)
        self.assertEqual(nnd[0][0], 6)
        self.assertEqual(nnd[0][1], 24)
        self.assertAlmostEqual(nnd[0][2], 1.6996583336037345)

    def test_nndist_nn(self):
        nnd = s.nearest_neighbors(self.cell, num_neigh=2)
        self.assertEqual(nnd[0][0], 0) 
        self.assertEqual(nnd[0][1], 29) 
        self.assertEqual(nnd[1][0], 0) 
        self.assertEqual(nnd[2][0], 1) 

    def test_coord_trans(self):
        self.assertEqual(array_almost_equal(
            np.linalg.solve(self.cell.get_cell().T, \
                            self.cell.get_positions().T).T % 1.0, \
            self.cell.get_scaled_positions()), True)
        

class RotateMolTestCase(unittest.TestCase):
    """Test the rotate_molecule function."""

    def test_rotate(self):
        coords = np.random.randn(10, 3)
        rcoords = s.rotate_molecule(coords, phi=math.pi, theta=2*math.pi, \
                psi=math.pi)
        self.assertEqual(array_almost_equal(coords, rcoords, 1e-15), True)

class InterpolateTestCase(unittest.TestCase):
    """Test the interpolate function."""

    def setUp(self):
        (self.c1, self.c2 ) = load_cells()

    def test_interpolate(self):
        c = s.interpolate_cells(self.c1, self.c2, frac=0.5)
        self.assertAlmostEqual(c.get_scaled_positions()[29,2], 0.258209253419)

    def test_interpolate_2(self):
        cells = s.interpolate_cells(self.c1, self.c2, images=2)
        self.assertAlmostEqual(cells[0].get_scaled_positions()[29,2], 0.24154259)
        self.assertAlmostEqual(cells[1].get_scaled_positions()[29,2], 0.27487592)

class AtomsMovedTestCase(unittest.TestCase):
    """Test the atoms_moved function."""

    def setUp(self):
        (self.c1, self.c2 ) = load_cells()

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
        (self.c1, self.c2 ) = load_cells()

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
    interpolate_suite = unittest.TestLoader().loadTestsFromTestCase(\
            InterpolateTestCase)
    return unittest.TestSuite([cell_suite, rot_suite, move_suite, \
            check_suite, interpolate_suite])


if __name__ == "__main__":
    unittest.main()
