#!/usr/bin/python

from distutils.core import setup

setup(name='vasputil',
      version='4.2',
      description='VASP utilities',
      author='Janne Blomqvist',
      author_email='Janne.Blomqvist@tkk.fi',
      url='http://www.fyslab.hut.fi/~/job/',
      packages=['vasputil'],
      package_dir = {'': 'lib/python'},
      scripts = ['bin/dosplot_example']
     )

