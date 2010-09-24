# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8
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

"""This module defines a class that represents a plane in 3d space."""

import numpy as np


class Plane(object):
    """Class for representing a plane in 3D space."""
    
    def __init__(self, points, normal=None):
        """Initialize plane.

        Arguments:
        points: points in the plane.
        normal: Normal vector of the plane.

        If normal is not provided, points must be a sequence or Numpy ndarray
        providing coordinates of 3 points in the plane. If coordinates are
        provided as a Numpy ndarray, each coordinate must be a row vector. If
        normal is provided, a single point is sufficient.

        """
        if normal == None:
            if isinstance(points, np.ndarray) and points.shape != (3,3):
                raise TypeError("Shape of points array must be (3,3).")
            elif len(points) != 3:
                raise TypeError("Points sequence must have 3 elemnents.")
            v1 = points[1] - points[0]
            v2 = points[2] - points[1]
            self.normal = np.cross(v1, v2)
            self.normal /= np.linalg.norm(self.normal)
            self.d_origo = np.dot(self.normal, points[1])
        else:
            self.normal = normal / np.linalg.norm(normal)
            self.d_origo = np.dot(self.normal, points)

    def distance(self, point):
        """Measure the distance between the plane and a point."""
        return abs(np.dot(self.normal, point) - self.d_origo)

# End of class Plane

def vec_pbc(vec):
    """Vector taking into account periodic boundary conditions.

    Input vector must be in direct coordinates.

    """
    arr = np.array(vec)
    v1 = arr.reshape(arr.size)
    v1 -= v1.round()
    return arr

def norm_pbc(vec):
    """Norm of a vector, taking into account periodic boundary conditions.

    Input vector must be in direct coordinates. If the input is a 2D array, it
    is assumed to be a list of vectors, and hence return an 1D array of the
    same length as the number of rows in the input array.

    """
    arr = np.array(vec)
    if len(arr.shape) == 1:
        return np.linalg.norm(arr - arr.round())
    elif len(arr.shape) == 2:
        nl = np.empty(arr.shape[0])
        for v in range(arr.shape[0]):
            nl[v] = np.linalg.norm(arr[v] - arr[v].round())
        return nl
    else:
        raise TypeError("Invalid shape of input")

