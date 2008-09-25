#!/usr/bin/python

from distutils.core import setup

setup(name='vasputil',
      version='4.3',
      description='VASP utilities',
      author='Janne Blomqvist',
      author_email='Janne.Blomqvist@tkk.fi',
      url='http://www.fyslab.hut.fi/~job/',
      packages=['vasputil'],
      package_dir = {'': 'lib/python'},
      scripts = ['bin/vasputil_atomsdistance', \
              'bin/vasputil_atomsmoved', \
              'bin/vasputil_direct2cartesian', \
              'bin/vasputil_dosplot', \
              'bin/vasputil_interpolate', \
              'bin/vasputil_nearestneighbors', \
              'bin/vasputil_poscar2xyz', \
              'bin/vasputil_zlayers']
     )

