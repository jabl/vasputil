# Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I config

MYSRC = src
SUBDIRS = $(MYSRC) bin

matlabdir = $(libdir)/m/
dist_matlab_DATA = $(srcdir)/lib/m/printcoords.m $(srcdir)/lib/m/rotatez.m

MYDOC = $(srcdir)/doc
robodocdir = $(prefix)/doc
dist_robodoc_DATA = $(MYDOC)/vasputil.html $(MYDOC)/vasputil.css

doc: all
	robodoc --src $(MYSRC) --doc $(MYDOC)/$(PACKAGE) --singledoc --html --sections --toc --nodesc

depend:
	cd $(MYSRC) && $(MAKE) $(AM_MAKEFLAGS) depend;
