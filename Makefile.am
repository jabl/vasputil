# Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I config

SUBDIRS = src

matlabdir = $(libdir)/m/
dist_matlab_DATA = $(srcdir)/lib/m/printcoords.m $(srcdir)/lib/m/rotate_molecule.m \
                   $(srcdir)/lib/m/translate_molecule.m

MYDOC = doc
robodocdir = $(prefix)/doc
dist_robodoc_DATA = $(MYDOC)/vasputil.html $(MYDOC)/vasputil.css

EXTRA_DIST = $(srcdir)/setup.py $(srcdir)/lib/python  $(srcdir)/bin

doc: $(MYDOC)/vasputil.html $(MYDOC)/vasputil.css

$(MYDOC)/vasputil.html $(MYDOC)/vasputil.css:
	mkdir doc
	robodoc --src $(srcdir)/src --doc $(MYDOC)/$(PACKAGE) --singledoc --html --sections --toc --nodesc

depend:
	cd src && $(MAKE) $(AM_MAKEFLAGS) depend;

install-data-local:
	@PYTHON@ ./setup.py install --install-purelib=$(pythondir) --install-scripts=$(bindir) 
