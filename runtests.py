#!/usr/bin/python
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
# Copyright (c) 2008 Janne Blomqvist

#  This file is part of Vasputil; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#  See the file COPYING for details.

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
