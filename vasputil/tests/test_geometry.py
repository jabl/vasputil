# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8

# Copyright (c) 2008, 2010 Janne Blomqvist

# This source code file is subject to the terms of the MIT (Expat)
# License. See the file LICENSE for details.

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
