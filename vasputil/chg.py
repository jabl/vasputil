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

"""This module defines a class that represents charge density. 

"""

try:
    import numpy as n
except ImportError:
    import pylab as n

import ase.io.vasp as aiv


class ChargeDensity(object):
    """Class for representing charge density."""
    
    def __init__(self, chg=None):
        """Initialize the data members of this class"""
        # VASP CHG/CHGCAR files contain supercell info as well.
        self.cell = None
        if chg != None:
            self.read_chg(chg)

    def read_chg(self, filename):
        """Read VASP CHG/CHGCAR file."""
        f = open(filename)
        self.atoms = aiv.read_vasp(f)
        f.readline()
        ng = f.readline().split()
        ng = (int(ng[0]), int(ng[1]), int(ng[2]))
        self.chg = n.empty(ng)
        # VASP writes charge density as
        # WRITE(IU,FORM) (((C(NX,NY,NZ),NX=1,NGXC),NY=1,NGYZ),NZ=1,NGZC)
        # Fortran nested implied do loops; innermost index fastest
        # First, just read it in
        for zz in range(self.chg.shape[2]):
            for yy in range(self.chg.shape[1]):
                self.chg[:, yy, zz] = n.fromfile(f, count = \
                        self.chg.shape[0], sep=' ')
        f.close()
        self.chg /= self.atoms.get_volume()

