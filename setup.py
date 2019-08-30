#!/usr/bin/env python3

from distutils.core import setup
from glob import glob

# Scripts whose names end in a-z or 1-9 (avoids emacs backup files)
scripts = glob('scripts/*[a-z,1-9]')

setup(name='vasputil',
      version='master',
      description='VASP utilities',
      author='Janne Blomqvist',
      author_email='Janne.Blomqvist@aalto.fi',
      url='https://github.com/jabl/vasputil',
      packages=['vasputil', 'vasputil.tests'],
      scripts = scripts)

