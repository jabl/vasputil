FC := ifort

OFLAG := -O2
DEBUG := -g -C
LIBS := -llapack -lf77blas -latlas
FFLAGS := $(DEBUG) -warn all -fpscomp none -std95
MODDIR := -module $(SRCDIR) -i $(SRCDIR)
