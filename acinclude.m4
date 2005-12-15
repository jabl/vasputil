m4_include([config/acx_f90.m4])
m4_include([config/acx_f90_libdir.m4])
m4_include([config/acx_f90_module.m4])
m4_include([config/acx_f90_modname.m4])
m4_include([config/acx_blas.m4])
m4_include([config/acx_lapack.m4])
m4_include([config/check_gnu_make.m4])
m4_include([config/acx_mpi.m4])
m4_include([config/acx_getkw.m4])

AC_DEFUN([ACX_BUILD_FLAGS],[. ./config/$1.conf])

AC_DEFUN([ACX_SUBST_BUILD_FLAGS],[
AC_SUBST(fopt, $fopt)
AC_SUBST(fwarn, $fwarn)
AC_SUBST(fdebug, $fdebug)
AC_SUBST(fprof, $fprof)
AC_SUBST(frange, $frange)
AC_SUBST(finclude, $finclude)
AC_SUBST(fldflags, $fldflags)
AC_SUBST(flibs, $flibs)
AC_SUBST(fdefs, $fdefs)

AC_SUBST(copt, $copt)
AC_SUBST(cwarn, $cwarn)
AC_SUBST(cdebug, $cdebug)
AC_SUBST(cprof, $cprof)
AC_SUBST(crange, $crange)
AC_SUBST(cinclude, $cinclude)
AC_SUBST(cldflags, $cldflags)
AC_SUBST(clibs, $clibs)
AC_SUBST(cdefs, $cdefs)

]) #ACX_SUBST_BUILD_FLAGS
