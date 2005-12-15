SYMLINKS = atomsdistance atomsmoved interpolate normalize unnormalize \
	planetoatom planetolayer POSCAR2xyz removeatoms \
	kspace2xyz check_nndist scgenerator lockatoms

bin_SCRIPTS = interpolate4 layers $(SYMLINKS)

bin_PROGRAMS = vasputil
vasputil_SOURCES = vasputil.f90 supercell_io.f90 supercell_core.f90 \
	geometry.f90 lattice.f90 mathconst.f90 kind_params.f90 conf.f90 \
	f2kcli.f90
	la_auxmod.f90 la_dgesdd.f90 la_dgesv.f90 \
	f95_lapack.f90 la_erinfo.f90 f77_lapack.f90

$(SYMLINKS):
	for syml in $(SYMLINKS); do\
		ln -sf vasputil $$syml;
	done

