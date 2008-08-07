#
# Process this file with autoconf to produce a configure script.
#
# Id: configure.ac 1.7 05/05/05 22:13:58+02:00 jonas@radon 
# $Id$

AC_INIT(vasputil, 4.1, Janne.Blomqvist@tkk.fi)
AC_CONFIG_SRCDIR(src/vasputil.F90)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE

AM_CONFIG_HEADER(src/config.h)

AC_PREFIX_DEFAULT([$HOME])
AC_PROG_MAKE_SET
AC_CANONICAL_HOST

# Check for programs.

CHECK_GNU_MAKE()
if test "$ifGNUmake" = "#"; then
	AC_MSG_ERROR([No working version of GNU make found!])
fi
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CHECK_PROGS(LD,ld,false)
AC_CHECK_PROGS(AR,ar,false)
AC_CHECK_PROGS(STRIP,strip,true)
AC_CHECK_PROGS(ETAGS,etags,true)
AC_CHECK_PROGS(CTAGS,ctags,true)
AC_PROG_LN_S

AC_CACHE_SAVE

AC_PROG_CC

# Trick from http://sources.redhat.com/ml/automake/2004-12/msg00070.html
#user_FCFLAGS = $FCFLAGS
AC_PROG_FC(pathf95 pathf90 pgf95 pgf90 ifort ifc lf95 xlf95 xlf90 gfortran f95 f90 fort)
#if test "x$user_FCFLAGS" = x; then
  # If the user didn't specify FCFLAGS, then FCFLAGS contains
  # a subset of -g -O2 selected by AC_PROG_FC.  This is not
  # a user setting, and we want to be able to override this
  # locally in our rules, so put these flags in a separate
  # variable and empty CFLAGS.
#  AC_SUBST([DEFAULTFLAGS], [$FCFLAGS])
#  FCFLAGS=
#fi


AC_ARG_ENABLE(debug-flags,
[AC_HELP_STRING([--enable-debug-flags@<:@=ARG@:>@],[compile with debugging flags (default is on)])
AC_HELP_STRING([--disable-debug-flags],[do not compile with debugging flags])
],USE_DEBUG_FLAGS="${enableval}",USE_DEBUG_FLAGS="yes")
AC_MSG_CHECKING([for debugging flags to use])

ACX_SUBST_FCFLAGS
AC_CACHE_SAVE

AX_FC_SEARCH_PATH

AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

#AC_PROG_YACC
#AC_PROG_LEX
#test "$YACC" = "bison -y" && YACC="$YACC -d"

dnl ACX_GETKW([], AC_MSG_ERROR([libgetkw.a not found!]))

AC_LANG_PUSH([Fortran])
ACX_BLAS(,AC_MSG_ERROR([No working blas library found!]))
ACX_LAPACK(,AC_MSG_ERROR([No working lapack library found!]))

AC_FC_SRCEXT([f90])
AC_FC_SRCEXT([F90])
#AC_FC_FREEFORM

#set -x -v
AX_F90_MODULE_EXTENSION
if test "x$ax_f90_modext" = "xunknown" ; then
     AC_MSG_ERROR([unable to find f90 modules extension])
fi
#set +x +v

AC_SUBST(AX_F90_MODEXT, $ax_f90_modext)

#AX_F90_MODULE_FLAG

# Check whether F2003 command line intrinsics are available
ACX_FC_COMMAND_LINE_ARGUMENTS

AC_LANG_POP([Fortran])

AC_CACHE_SAVE
#AC_MSG_NOTICE([Checking if MPI is available])
#AC_LANG_PUSH([Fortran 90])
#ACX_MPI([FC=$MPIFC])
#AC_LANG_POP
#fdefs="$fdefs -DHAVE_MPI"
#AC_CACHE_SAVE

#AC_SUBST(top_bindir,['$(topdir)/bin'])
#AC_SUBST(top_includedir,['$(topdir)/include'])
#AC_SUBST(top_libdir,['$(topdir)/lib'])

AC_CONFIG_FILES([Makefile src/Makefile bin/Makefile])
AC_OUTPUT

ACX_F90_CONFIG_HEADER([src/config.h], [src/config.inc])
