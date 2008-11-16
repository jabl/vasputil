#!/usr/bin/python

from distutils.core import setup
from glob import glob

scripts = glob('scripts/*')

setup(name='vasputil',
      version='5.0',
      description='VASP utilities',
      author='Janne Blomqvist',
      author_email='Janne.Blomqvist@tkk.fi',
      url='http://www.fyslab.hut.fi/~job/',
      packages=['vasputil'],
      scripts = scripts)

