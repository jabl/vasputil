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

class CellTestCase(unittest.TestCase):
    """Testcase for vasputil.supercell.Cell class."""

    def setUp(self):
        # WARNING: Ugly ugly ugly!
        # find the path where this file is
        path = os.path.split(__file__)[0]
        # The test POSCAR file is in the same directory.
        path = os.path.join(path, "POSCAR")
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
        dist = self.cell.atoms_distance(9, 24, "xy")
        ntu.assert_almost_equal(dist, 1.51230705052) 

class RotateMolTestCase(nt.NumpyTestCase):
    """Test the rotate_molecule function."""

    def test_rotate(self):
        coords = ntu.rand(10, 3)
        rcoords = s.rotate_molecule(coords, phi=math.pi, theta=2*math.pi, \
                psi=math.pi)
        ntu.assert_array_almost_equal(coords, rcoords)

def suite():
    cell_suite = unittest.TestLoader().loadTestsFromTestCase(CellTestCase)
    rot_suite = unittest.TestLoader().loadTestsFromTestCase(RotateMolTestCase)
    return unittest.TestSuite([cell_suite, rot_suite])


if __name__ == "__main__":
    unittest.main()
