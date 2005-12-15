dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl 	$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F90_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F90 libraries.
dnl Users will also need to use AC_F90_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F90 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
dnl Modified by Jonas Juselius <jonas@iki.fi>
dnl
AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])

acx_lapack_save_LIBS="$LIBS"
acx_lapack_save_LDFLAGS="$LDFLAGS"
acx_lapack_libs=""
acx_lapack_dir=""

AC_ARG_WITH(lapack,
	[AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])

case $with_lapack in
	yes | "") ;;
	no) acx_lapack_ok=disable ;;
	-l* | */* | *.a | *.so | *.so.* | *.o) acx_lapack_libs="$with_lapack" ;;
	*) acx_lapack_libs="-l$with_lapack" ;;
esac

AC_ARG_WITH(lapackdir,
	[AC_HELP_STRING([--with-lapackdir=<dir>], [look for LAPACK library in <dir>])])

case $with_lapackdir in
     "") ;;
	-L*) LDFLAGS="$LDFLAGS $with_lapackdir" 
	     acx_lapack_dir="$with_lapackdir" ;;
	*) LDFLAGS="$LDFLAGS -L$with_lapackdir " 
	   acx_lapack_dir="-L$with_lapackdir" ;;
esac

if test "x$acx_lapack_ok" != xyes; then
	acx_lapack_ok=no
fi

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
	acx_lapack_ok=noblas
fi

# Get fortran linker name of LAPACK function to check for.
AC_F90_FUNC(cheev)

# If --with-lapack is defined, then look for THIS AND ONLY THIS lapack lib
if test $acx_lapack_ok = no; then
case $with_lapack in
    ""|yes) ;;
	*) save_LIBS="$LIBS"; LIBS="$acx_lapack_libs $BLAS_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $cheev in $acx_lapack_libs])
	AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=specific])
	AC_MSG_RESULT($acx_lapack_ok)
	LIBS="$save_LIBS"
	;;
esac
fi


# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $cheev in $LAPACK_LIBS])
	AC_TRY_LINK_FUNC($cheev, [acx_lapack_ok=yes; acx_lapack_libs=$LAPACK_LIBS])
	AC_MSG_RESULT($acx_lapack_ok)
	LIBS="$save_LIBS"
fi

# Intel mkl LAPACK library?
if test $acx_lapack_ok = no; then
	AC_CHECK_LIB(mkl_lapack, $cheev,
		[acx_lapack_ok=yes; acx_lapack_libs="-lmkl_lapack"], [], 
		[$BLAS_LIBS $FLIBS])
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
	AC_CHECK_FUNC($cheev, [acx_lapack_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic LAPACK library?
if test $acx_lapack_ok = no; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_CHECK_LIB(lapack, $cheev,
		[acx_lapack_ok=yes; acx_lapack_libs="-llapack"], [], [$FLIBS])
	LIBS="$save_LIBS"
fi

LAPACK_LIBS="$acx_lapack_libs" 
AC_SUBST(LAPACK_LIBS)

LIBS="$acx_lapack_save_LIBS"
LDFLAGS="$acx_lapack_save_LDFLAGS $acx_lapack_dir"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK
