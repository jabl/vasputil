#!/usr/bin/pythoon

# Plot DOS from vasp
# First use split_dos to create DOSx files.

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
