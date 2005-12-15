FC := pathf90

OFLAG := -Ofast
DEBUG := -g -O0 -C -fno-unsafe-math-optimizations -msse2 \
	-TENV:simd_omask=OFF -TENV:simd_umask=OFF -TENV:simd_imask=OFF \
	-TENV:simd_dmask=OFF -TENV:simd_zmask=OFF -TENV:simd_pmask=OFF
LIBS := /opt/acml/pathscale64/lib/libacml.a
FFLAGS := $(DEBUG) -ansi -fullwarn -Wall
MODDIR := -module $(SRCDIR) -I $(SRCDIR)
