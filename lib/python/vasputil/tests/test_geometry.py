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

"""This module contains unit tests for the vasputil.geometry module."""

import unittest
import vasputil.geometry as g
import numpy.testing.numpytest as nt
import numpy.testing.utils as ntu

class PlaneTestCase(nt.NumpyTestCase):
    """Testcase for vasputil.geometry.Plane class."""

    def test_construct(self):
        self.assertRaises(TypeError, g.Plane)

    def test_construct_normal(self):
        pl1 = g.Plane((127.0, 37.0, 42.0), (0., 0., 10.))
        ntu.assert_almost_equal(pl1.d_origo, 42.0)

def suite():
    plane_suite = unittest.TestLoader().loadTestsFromTestCase(PlaneTestCase)
    return unittest.TestSuite([plane_suite])


if __name__ == "__main__":
    unittest.main()
