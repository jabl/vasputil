# -*- coding: latin-1 -*-
# vim: set fileencoding=latin-1

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

"""This module contains unit tests for the vasputil.geometry module."""

import unittest
import vasputil.geometry as g

class PlaneTestCase(unittest.TestCase):
    """Testcase for vasputil.geometry.Plane class."""

    def test_construct(self):
        self.assertRaises(TypeError, g.Plane)

    def test_construct_normal(self):
        pl1 = g.Plane((127.0, 37.0, 42.0), (0., 0., 10.))
        self.assertAlmostEqual(pl1.d_origo, 42.0)

class NormPbcTestCase(unittest.TestCase):
    """Test the norm_pbc function."""

    def test_norm_pbc(self):
        a = [1.5, 0, 0]
        d = g.norm_pbc(a)
        self.assertAlmostEqual(d, 0.5)

    def test_norm_pbc_2(self):
        a = [1.0, 1.0, 1.0]
        d = g.norm_pbc(a)
        self.assertAlmostEqual(d, 0.0)

def suite():
    plane_suite = unittest.TestLoader().loadTestsFromTestCase(PlaneTestCase)
    norm_suite = unittest.TestLoader().loadTestsFromTestCase(NormPbcTestCase)
    return unittest.TestSuite([plane_suite, norm_suite])


if __name__ == "__main__":
    unittest.main()
