bin_SCRIPTS = layers
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = layers.in

do_subst = sed -e 's,[@]PACKAGE_STRING[@],$(PACKAGE_STRING),g'

layers: layers.in Makefile
	$(do_subst) < $(srcdir)/layers.in > layers
	chmod +x layers
