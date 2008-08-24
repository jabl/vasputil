#!/usr/bin/pythoon
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



"""Utility function to set some default plot parameters for plotting DOS
figures.

"""

from pylab import *


def showdp(fsz=16):
    xlabel("E-E$_\mathrm{f}$ (eV)", size=fsz)
    figtext(0.03, 0.45, "LDOS", rotation='vertical', size=fsz)
    loc, lab = xticks()
    lab.set_size = fsz
    loc, lab = yticks()
    lab.set_size = fsz
    legend()
    subplots_adjust(hspace=0.0)
    show()
