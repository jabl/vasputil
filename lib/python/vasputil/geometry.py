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

"""This module defines a class that represents a plane in 3d space."""

try:
    import numpy as n
except:
    import pylab as n


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
            if type(points) == n.ndarray and points.shape != (3,3):
                raise TypeError("Shape of points array must be (3,3).")
            elif len(points) != 3:
                raise TypeError("Points sequence must have 3 elemnents.")
            v1 = points[1] - points[0]
            v2 = points[2] - points[1]
            self.normal = n.cross(v1, v2)
            self.normal /= n.linalg.norm(self.normal)
            self.d_origo = n.dot(self.normal, points[1])
        else:
            self.normal = normal / n.linalg.norm(normal)
            self.d_origo = n.dot(self.normal, points)
