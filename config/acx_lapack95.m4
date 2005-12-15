dnl
dnl @synopsis ACX_LAPACK95([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Written by Jonas Juselius <jonas@iki.fi>
dnl

AC_DEFUN([ACX_LAPACK95], [
AC_PREREQ(2.57)
AC_REQUIRE([AC_PROG_F90])
AC_REQUIRE([ACX_LAPACK])

acx_lapack95_ok=no
acx_lapack95_save_FFLAGS="$FFLAGS"
acx_lapack95_save_LIBS="$LIBS"
acx_lapack95_save_LDFLAGS="$LDFLAGS"
acx_lapack95_paths="/usr/lib/lapack95 /usr/local/lib/lapack95 /opt/lib/lapack95 ./lapack95"

acx_lapack95_libs=""
acx_lapack95_dir=""
acx_lapack95_includes=""

AC_ARG_WITH([lapack95], AC_HELP_STRING([--with-lapack95=<lib>]))

AC_ARG_WITH([lapack95dir], AC_HELP_STRING([--with-lapack95dir=/usr/local/lib],
	[base path to the lapack95 installation]))

case $with_lapack95dir in
	  yes | no | "") ;;
	  *) LDFLAGS="$LDFLAGS -L$with_lapack95dir" 
		 FFLAGS="$FFLAGS -I$with_lapack95dir"
		 acx_lapack95_dir="-L$with_lapack95dir" 
		 acx_lapack95_includes="-I$with_lapack95dir" ;;
esac

case $with_lapack95 in
	yes | "") ;;
	no) acx_lapack95_ok=disable ;;
	-l* | */* | *.a | *.so | *.so.* | *.o) 
		acx_lapack95_libs="$with_lapack95" ;;
	*) acx_lapack95_libs="-l$with_lapack95" ;;
esac

acx_lapack95_mod_ok=no

#
# Test if f95_lapack.mod is available
#
AC_LANG_PUSH([Fortran 90])
AC_MSG_CHECKING([for f95_lapack.mod])
AC_COMPILE_IFELSE(AC_LANG_PROGRAM([],
	[use f95_lapack]), 
   [acx_lapack95_mod_ok=yes], [acx_lapack95_includes=""])
if test $acx_lapack95_mod_ok = no; then
save_FFLAGS=$FFLAGS
for i in $acx_lapack95_paths; do
  FFLAGS="-I$i $save_FFLAGS"
  AC_COMPILE_IFELSE(AC_LANG_PROGRAM([],[use f95_lapack]), 
	[acx_lapack95_mod_ok=yes; acx_lapack95_includes=-I$i])
   if test $acx_lapack95_mod_ok = yes; then
	break
   fi
done
FFLAGS=$save_FFLAGS
fi
AC_MSG_RESULT($acx_lapack95_mod_ok)

# Get fortran linker names of LAPACK95 functions to check for.
AC_F90_FUNC(sgeev_f95)

#
# Test if lapack95 is linkable
#
# If --with-lapack95 is defined, then look for THIS AND ONLY THIS lapack lib
if test "x$with_lapack95" != x; then
if test $acx_lapack95_mod_ok = yes; then
	save_LIBS="$LIBS"
	LIBS="$acx_lapack_libs $LAPACK_LIBS $BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgeev_f95 in $acx_lapack95_libs])
	AC_TRY_LINK_FUNC(sgeev_f95, [acx_lapack95_ok=yes], [acx_lapack95_ok=specific])
	AC_MSG_RESULT($acx_lapack95_ok)
	LIBS="$save_LIBS"
fi; fi

if test $acx_lapack95_mod_ok = yes; then
save_LIBS="$LIBS"
LIBS="-llapack95 $LAPACK_LIBS $BLAS_LIBS $LIBS"
AC_MSG_CHECKING([for $sgeev_f95 in liblapack95])
AC_TRY_LINK_FUNC(sgeev_f95, [acx_lapack95_ok=yes;
    acx_lapack95_libs=-llapack95])
if test $acx_lapack95_ok = no; then
	for i in $acx_lapack95_paths; do
	  save_LDFLAGS=$LDFLAGS
	  LDFLAGS="-L$i $save_LDFLAGS"
	  AC_TRY_LINK_FUNC(sgeev_f95, [acx_lapack95_ok=yes;
	  acx_lapack95_dir=-L$i; acx_lapack95_libs=-llapack95])
	  if test $acx_lapack95_ok = yes; then
		break
	  fi
	  LDFLAGS=$save_LDFLAGS
	done
	LIBS="$save_LIBS"
fi
AC_MSG_RESULT($acx_lapack95_ok)
fi
AC_LANG_POP()

LIBS="$acx_lapack95_save_LIBS"
FFLAGS="$acx_lapack95_includes $acx_lapack95_save_FFLAGS"
LDFLAGS="$acx_lapack95_dir $acx_lapack95_save_LDFLAGS"

LAPACK_LIBS="$acx_lapack95_libs $LAPACK_LIBS"
AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test $acx_lapack95_ok = no; then
		if test "x$acx_lapack95_ok" = "xdisable"; then
			acx_lapack95_ok=no
		else
			acx_lapack95_ok=no
			$2
		fi
else
		ifelse([$1],,AC_DEFINE(HAVE_LIBLAPACK95,1,
		[Define if you have the lapack95 library.]),[$1])
		:
fi
]) dnl ACX_LAPACK95
