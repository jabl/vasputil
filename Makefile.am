# Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I config

MYSRC = $(srcdir)/src
SUBDIRS = src bin

matlabdir = $(libdir)/m/
dist_matlab_DATA = $(srcdir)/lib/m/printcoords.m $(srcdir)/lib/m/rotate_molecule.m \
                   $(srcdir)/lib/m/translate_molecule.m

MYDOC = doc
robodocdir = $(prefix)/doc
dist_robodoc_DATA = $(MYDOC)/vasputil.html $(MYDOC)/vasputil.css

#mypydir = $(srcdir)/lib/python/vasputil
#python_PYTHON = $(srcdir)/setup.py $(mypydir)/__init__.py $(mypydir)/dosplot.py \
#		   $(mypydir)/supercell.py

EXTRA_DIST = $(srcdir)/lib/python

doc: $(MYDOC)/vasputil.html $(MYDOC)/vasputil.css

$(MYDOC)/vasputil.html $(MYDOC)/vasputil.css:
	mkdir doc
	robodoc --src $(MYSRC) --doc $(MYDOC)/$(PACKAGE) --singledoc --html --sections --toc --nodesc

depend:
	cd src && $(MAKE) $(AM_MAKEFLAGS) depend;

install-data-local:
	@PYTHON@ ./setup.py install --install-purelib=$(pythondir) --install-scripts=$(bindir) 
