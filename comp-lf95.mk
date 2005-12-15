FC := lf95

OFLAG := --o2 --ntrace --nap --nchk --ng --nsav --prefetch 2 --tpp
DEBUG := -g --ap --chkglobal --lst --trace --pca --xref
LIBS := -llapackmt -lblasmt
LDFLAGS := --staticlink
FFLAGS := $(DEBUG) --f95 --info --wo --warn --trap
MODDIR := -M $(SRCDIR)
