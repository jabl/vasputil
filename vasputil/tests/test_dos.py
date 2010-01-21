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

"""This module contains unit tests for the vasputil.dos module."""

import unittest
import vasputil.dos as d

class LdosTestCase(unittest.TestCase):
    """Testcase for vasputil.dos.LDOS class."""


def suite():
    ldos_suite = unittest.TestLoader().loadTestsFromTestCase(LdosTestCase)
    return unittest.TestSuite([ldos_suite])


if __name__ == "__main__":
    unittest.main()
