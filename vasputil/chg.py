# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
# Copyright (c) 2010 Janne Blomqvist

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

"""Some helper functions for handling charge densities.
"""

import numpy as np

def load_chg_plane(fname):
    """Load charge density in a plane, in the format used by lev00"""
    dat = np.loadtxt(fname)

    # Data is stored in text form, with each line containing
    # X Y Z
    # where Z is the charge density
    # However, if the width along the axises are different, the grid points in
    # each direction stay constant, but the grid spacing changes.
    # So we might need to reinterpolate on a regular grid

    # Original num of points in both X and Y directions.
    n = np.sqrt(dat.shape[0]) 
    # x, y = meshgrid(dat[0:n, 0], dat[0:n, 1])
    x = dat[:,0]
    y = dat[:,1]
    xmax = x.max()
    xmin = x.min()
    ymax = y.max()
    ymin = y.min()
    xr = xmax - xmin
    yr = ymax - ymin

    xl = np.linspace(xmin, xmax, n)
    yl = np.linspace(ymin, ymax, n)
    z = np.reshape(dat[:, 2], (n, n))

    if xr == yr:
        # No interpolation needed
        pass
    else:
        import scipy.interpolate as spi
        # Create 2D interpolator object with original data
        # interp2d goes into an infinite loop?
        # inter2 = spi.interp2d(xl, yl, z)
        inter2 = spi.RectBivariateSpline(xl, yl, z)
        # Reinterpolate only along the bigger dimension
        idim = xr < yr
        if idim: # Y range is bigger
            fac = int(round(n * yr / xr))
            yl = np.linspace(ymin, ymax, fac)
        else:
            fac = int(round(n * xr / yr))
            xl = np.linspace(xmin, xmax, fac)
        z = inter2(xl, yl)

    return xl, yl, z
