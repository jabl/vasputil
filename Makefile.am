# Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS = -I config

SUBDIRS = src

MYDOC = doc
apidoc = fortrandoc
docdir = ${datadir}/doc/${PACKAGE}

dist_doc_DATA = $(MYDOC)/$(apidoc).html $(MYDOC)/$(apidoc).css $(MYDOC)/README.html

EXTRA_DIST = $(srcdir)/setup.py $(srcdir)/lib/python  $(srcdir)/bin

doc: $(MYDOC)/$(apidoc).html $(MYDOC)/$(apidoc).css $(MYDOC)/README.html

$(MYDOC)/$(apidoc).html $(MYDOC)/$(apidoc).css:
	mkdir -p doc
	robodoc --src $(srcdir)/src --doc $(MYDOC)/$(apidoc) --singledoc --html --sections --toc --nodesc

$(MYDOC)/README.html: $(srcdir)/README
	mkdir -p doc
	rst2html $(srcdir)/README $(MYDOC)/README.html

depend:
	cd src && $(MAKE) $(AM_MAKEFLAGS) depend;

test:
	@PYTHON@ $(srcdir)/lib/python/runtests.py

install-data-local:
	@PYTHON@ ./setup.py install --install-purelib=$(pythondir) --install-scripts=$(bindir) 
