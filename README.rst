========
Vasputil
========

Vasputil is a set of command-line utilities and Python libraries
designed to make life with VASP easier. The command-line utilities
(Python scripts), can be used directly or as examples of what can be
done with the provided python modules.  These are

``vasputil_atomsdistance``
    Measure distance between two atoms, can also measure projected
    distances.

``vasputil_atomsmoved`` 
    See which atoms have moved between two POSCAR files.  If lattices are
    compatible, takes into account periodic boundary conditions.  Optionally
    print out moved distances.

``vasputil_chgarith``
    Simple arithmetic on charge density in CHG/CHGCAR format
    files. Supports +,-,*,/ for the elemental arithmetic operations
    and 'avg' for calculating the average.

``vasputil_chgplaneplot``
    Plot charge density in a plane (generated e.g. by lev00) as a
    pseudocolor, contour (colour or BW), or filled contour
    plot. Output is either to a window on the screen or to a file.

``vasputil_direct2cartesian``
    Convert a POSCAR file from direct to cartesian coordinates.

``vasputil_dosplot``
    Plot density-of-states. Can plot total, total integrated, and
    site-projected DOS. For site-projected DOS, different orbitals can
    be plotted.

``vasputil_interpolate``
    Interpolate coordinates between two POSCAR files, either one intermediate
    image where one can specify where in the interval the new coordinates will
    be, or create many evenly spaced images.

``vasputil_nearestneighbors``
    Print out nearest neighbor table, based on minimum distance or number of
    nearest neighbors.

``vasputil_plane2atom``
    Calculate the distance between an atom and a plane defined by three atoms.

``vasputil_poscar2xyz``
    Convert a POSCAR format file to XYZ format.

``vasputil_zlayers``
    Find layers in the z direction and print interlayer distances.

Documentation for these utilities is provided via the ``-h`` command
line option.

Web Page
--------

http://github.com/jabl/vasputil
    Main home page; download, source tree, issue tracker, etc.

http://tfy.tkk.fi/~job/vasputil/
    Older releases, and a copy of the README page.

Requirements
------------

For full functionality, vasputil requires the following python libraries:

1) `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`_ (Atomic
   Simulation Environment) version 3.4. For older ASE releases,
   vasputil 5.3 is compatible with ASE 3.3.1, vasputil 5.2 is
   compatible with ASE 3.2, vasputil 5.1 with ASE 3.1, and vasputil
   5.0 with ASE 3.0.

2) `NumPy <http://www.scipy.org/NumPy>`_ 

3) `matplotlib <http://matplotlib.sf.net>`_

4) `SciPy <http://www.scipy.org/SciPy>`_

Finally, to regenerate the README.html file, the `rst2html
<http://docutils.sourceforge.net/docs/user/tools.html>`_ utility is
needed.

Without all the libraries, a subset of the functionality is
available. The matplotlib library is needed only by the
vasputil_chgplaneplot, vasputil_dosplot and vasputil_zlayers
utilities. The SciPy library is needed only by vasputil_chgplaneplot,
and only in case the figure shape is not square. The remaining
utilities require only ASE and NumPy.

Except for ASE, the other three libraries (and their dependencies)
should already be installed on your system. If not, at least on Linux
they should be available from the default package repositories, and
installation should be a breeze.  If you don't have access to a
package management system (e.g. on Linux) which includes matplotlib
and its dependencies, you can get all at once by installing one of the
python distributions targeted at science, such as SAGE, Enthought
Python Distribution, or Python(x,y). Another solution if you already
have python installed, is to install matplotlib and other libraries
via *easy_install* from setuptools (a sort of package management
system for python packages).


Installation Instructions
-------------------------

1)  Unpack the program distribution and go to the directory.

    ``tar xvzf vasputil-$(VERSION).tar.gz``

    ``cd vasputil-$(VERSION)``

2)  Install the program into your home directory

    ``python setup.py install --home=~``

    This will install the scripts into ~/bin and the python modules
    into ~/lib/python/vasputil. You should make sure that the
    directory of the python modules is found on the python search
    path. E.g. if you installed it like the above, you should have

    ``export PYTHONPATH=$HOME/lib/python``

    in your shell initialization script (~/.bash_profile or
    equivalent). For installation into some other directory, see
    ``python setup.py --help`` and the documentation for `distutils
    <http://www.python.org/doc/lib/module-distutils.html>`_ .

    For SUSE Linux users: There is apparently something weird with the
    python installation in SUSE, and the ``--home`` option does not
    work properly. Instead install vasputil with

    ``python setup.py install --prefix=~``

    In this case the library modules will instead be installed in
    ~/lib{64}/python{X.Y}/site-packages/vasputil.

3)  Optionally, run the testsuite to make sure it works.

    ``python runtests.py``

Interactive Usage
-----------------

The python modules are designed to be used from the interactive python prompt
as well, to provide a sort of MATLAB-like interactive environment where you can
manipulate coordinates stored in arrays, or plot data etc. For this, the
`IPython <http://ipython.scipy.org/moin/>`_ enhanced interactive environment is
recommended. Start with the matplotlib stuff preloaded with ``ipython -pylab``.
For those familiar with MATLAB, see `NumPy for MATLAB Users
<http://www.scipy.org/NumPy_for_Matlab_Users>`_ . After that, just import what
you need from the vasputil module and get going! See the python scripts for
reference, or ``pydoc vasputil``, ``pydoc vasputil.supercell`` etc. to read
documentation, or the ``?`` operator in IPython.

With IPython you can also run the scripts in the interactive
environment with the ``%run`` command, and this will import the
variables in the script into the ipython environment. This is useful
if you want to do something similar, but not quite what the
command-line interface allows.


Manual
------

This section describes how to accomplish specific tasks. It does not
document every feature of vasputil, as most of the tasks are hopefully
simple enough that the usage should be self-evident from the help
instruction given by the ``-h`` option to the command-line utilities.

Plotting charge density in a plane
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the ``vasputil_chgarith`` tool to create a suitable input charge
density file (e.g. the charge density difference when some species
adsorbs on a surface), and then use the `lev00
<http://www.cmmp.ucl.ac.uk/~lev/codes/lev00/>`_ utility to create the
charge density data in a plane. The data is written to a file (default
name out.dat_1) that can be read with the ``vasputil_chgplaneplot``
utility which can then plot this data.

Generating supercells
~~~~~~~~~~~~~~~~~~~~~

vasputil 5.x does not contain the supercell generator previously found
in vasputil 4.x. To create supercells, the ASE supercell generator is
recommended. This contains functionality to create surface supercells
as well. In order to generate a supercell for VASP, first create the
supercell following the `ASE lattice
<https://wiki.fysik.dtu.dk/ase/ase/lattice.html>`_ instructions. This
is probably easiest done interactively using IPython. Assuming you
have imported ase as ``from ase import *`` and your supercell is
referenced via the variable ``atoms``, you can visualize your results
with ``view(atoms)``. Finally, write out the supercell to the file
``POSCAR.out`` with ``write('POSCAR.out', atoms, format='vasp')``.

Another option is the `tetr
<http://www.cmmp.ucl.ac.uk/~lev/codes/lev00/>`_ utility, which also
contains a supercell generator.

Plotting projected Density-of-States
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``vasputil_dosplot`` contains a simple command-line utility for
plotting a single orbital. While this is nice for quickly looking at
the DOS, for publication plots you probably want to look at that
utility and create custom versions of it in order to create the
specific plots you want, with multiple subplots and multiple graphs
per plot.
