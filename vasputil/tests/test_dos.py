# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8

# Copyright (c) 2008, 2010 Janne Blomqvist

# This source code file is subject to the terms of the MIT (Expat)
# License. See the file LICENSE for details.

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
