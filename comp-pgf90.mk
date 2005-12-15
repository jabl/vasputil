FC := pgf90

OFLAG := -fast -Mcache_align
DEBUG := -g -Mbounds -Mchkfpstk 
LIBS := /opt/acml/pgi64/lib/libacml.a
FFLAGS := $(DEBUG) -Mstandard -Mchkptr 
MODDIR := -module $(SRCDIR)
