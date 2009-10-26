#!/usr/bin/python

from distutils.core import setup
from glob import glob

# Scripts whose names end in a-z or 1-9 (avoids emacs backup files)
scripts = glob('scripts/*[a-z,1-9]')

setup(name='vasputil',
      version='5.3',
      description='VASP utilities',
      author='Janne Blomqvist',
      author_email='Janne.Blomqvist@tkk.fi',
      url='http://tfy.tkk.fi/~job/vasputil',
      packages=['vasputil', 'vasputil.tests'],
      scripts = scripts)

