#!/usr/bin/python
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
# Copyright (c) 2008 Janne Blomqvist

# This source code file is subject to the terms of the MIT (Expat)
# License. See the file LICENSE for details.

"""This module imports all the vasputil unit tests and runs them all."""

import unittest
from vasputil.tests import *

def suite():
    geo_suite = test_geometry.suite()
    scell_suite = test_supercell.suite()
    dos_suite = test_dos.suite()
    return unittest.TestSuite([geo_suite, scell_suite, dos_suite])


if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(suite())
