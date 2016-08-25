# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
# Copyright (c) 2010, 2016 Janne Blomqvist

# This source code file is subject to the terms of the LGPL 2.1
# License. See the file LICENSE for details.

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
