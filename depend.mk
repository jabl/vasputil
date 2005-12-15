FOBJ=src/conf.o src/f2kcli.o src/f77_lapack.o src/f95_lapack.o src/geometry.o src/kind_params.o src/la_auxmod.o src/la_dgesdd.o src/la_dgesv.o src/la_erinfo.o src/la_precision.o src/lattice.o src/mathconst.o src/supercell_core.o src/supercell_io.o src/vasputil.o 

vasputil: $(FOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(FOBJ) $(LIBS)

src/conf.o : src/conf.f90 src/kind_params.o 
src/f2kcli.o : src/f2kcli.f90 
src/f77_lapack.o : src/f77_lapack.f90 src/la_precision.o 
src/f95_lapack.o : src/f95_lapack.f90 src/la_precision.o 
src/geometry.o : src/geometry.f90 src/f95_lapack.o src/conf.o 
src/kind_params.o : src/kind_params.f90 
src/la_auxmod.o : src/la_auxmod.f90 
src/la_dgesdd.o : src/la_dgesdd.f90 src/f77_lapack.o src/la_auxmod.o src/la_precision.o 
src/la_dgesv.o : src/la_dgesv.f90 src/f77_lapack.o src/la_auxmod.o src/la_precision.o 
src/la_erinfo.o : src/la_erinfo.f90 
src/la_precision.o : src/la_precision.f90 
src/lattice.o : src/lattice.f90 src/f95_lapack.o src/mathconst.o src/geometry.o src/conf.o 
src/mathconst.o : src/mathconst.f90 src/conf.o 
src/supercell_core.o : src/supercell_core.f90 src/lattice.o src/geometry.o src/conf.o 
src/supercell_io.o : src/supercell_io.f90 src/supercell_core.o src/conf.o 
src/vasputil.o : src/vasputil.f90 src/kind_params.o src/supercell_io.o src/f2kcli.o src/conf.o 
