bin_SCRIPTS = interpolate4 layers
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = interpolate4.in layers.in

do_subst = sed -e 's,[@]PACKAGE_STRING[@],$(PACKAGE_STRING),g'

interpolate4: interpolate4.in Makefile
	$(do_subst) < $(srcdir)/interpolate4.in > interpolate4
	chmod +x interpolate4

layers: layers.in Makefile
	$(do_subst) < $(srcdir)/layers.in > layers
	chmod +x layers
