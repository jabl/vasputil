# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8

# Copyright (c) 2008, 2010 Janne Blomqvist

#  This file is part of Vasputil; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#  See the file COPYING for details.

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
