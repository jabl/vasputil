conf.o : conf.f90 kind_params.o 
f2kcli.o : f2kcli.f90 
f77_lapack.o : f77_lapack.f90 la_precision.o 
f95_lapack.o : f95_lapack.f90 la_precision.o 
geometry.o : geometry.f90 f95_lapack.o conf.o 
kind_params.o : kind_params.f90 
la_auxmod.o : la_auxmod.f90 
la_dgesdd.o : la_dgesdd.f90 f77_lapack.o la_auxmod.o la_precision.o 
la_dgesv.o : la_dgesv.f90 f77_lapack.o la_auxmod.o la_precision.o 
la_erinfo.o : la_erinfo.f90 
la_precision.o : la_precision.f90 
lattice.o : lattice.f90 f95_lapack.o mathconst.o geometry.o conf.o 
mathconst.o : mathconst.f90 conf.o 
poscar_io.o : poscar_io.f90 supercell_utils.o supercell_core.o conf.o 
sc_file_convert.o : sc_file_convert.f90 supercell_utils.o poscar_io.o conf.o 
supercell_core.o : supercell_core.f90 lattice.o geometry.o conf.o 
supercell_generator.o : supercell_generator.f90 poscar_io.o supercell_utils.o supercell_core.o lattice.o geometry.o conf.o 
supercell_measure.o : supercell_measure.f90 supercell_utils.o poscar_io.o lattice.o geometry.o conf.o 
supercell_modify.o : supercell_modify.f90 supercell_utils.o poscar_io.o supercell_core.o conf.o 
supercell_utils.o : supercell_utils.f90 supercell_core.o lattice.o geometry.o conf.o 
vasputil.o : vasputil.f90 kind_params.o supercell_generator.o sc_file_convert.o supercell_modify.o supercell_measure.o f2kcli.o conf.o 
