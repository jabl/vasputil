FC := xlf95

OFLAG := -O3 -qstrict -qarch=pwr4 -qtune=pwr4 -qcache=auto
DEBUG := -g -C
FFLAGS := $(DEBUG) -qsuffix=f=f90 -qlanglvl=95pure
MODDIR := -qmoddir=$(SRCDIR)
LIBS := -lessl -llapack
