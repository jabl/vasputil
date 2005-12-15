dnl
dnl @synopsis ACX_GETKW([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Written by Jonas Juselius <jonas@iki.fi>
dnl

AC_DEFUN([ACX_GETKW], [
AC_PREREQ(2.57)
AC_REQUIRE([AC_PROG_F90])

acx_getkw_ok=no
acx_getkw_save_FFLAGS="$FFLAGS"
acx_getkw_save_CFLAGS="$CFLAGS"
acx_getkw_save_LIBS="$LIBS"
acx_getkw_save_LDFLAGS="$LDFLAGS"

acx_getkw_libs=""
acx_getkw_dir=""
acx_getkw_includes=""

AC_ARG_WITH([getkw], AC_HELP_STRING([--with-getkw=LIB],[]))

AC_ARG_WITH([getkw_includes], 
  AC_HELP_STRING([--with-getkw-modules=DIR],[]))
  
AC_ARG_WITH([getkw_libdir], 
  AC_HELP_STRING([--with-getkw-libdir=DIR],[]))

AC_ARG_WITH([getkwdir], AC_HELP_STRING([--with-getkwdir=DIR],
	[Base path to the libgetkw installation]))

case $with_getkwdir in
	  yes | no | "") ;;
	  *) acx_getkw_dir="-L$with_getkwdir/lib" 
		 acx_getkw_includes="-I$with_getkwdir/include" 
		 LDFLAGS="$LDFLAGS $acx_getkw_dir" 
		 FFLAGS="$FFLAGS $acx_getkw_includes"
		 CFLAGS="$CFLAGS $acx_getkw_includes"
	  ;;
esac

case $with_getkw_includes in
	yes | no | "") ;;
	*) acx_getkw_includes="-I$with_getkw_includes" 
	   FFLAGS="$FFLAGS $acx_getkw_includes"
	   CFLAGS="$CFLAGS $acx_getkw_includes"
	 ;;
esac

case $with_getkw_libdir in
	yes | no | "") ;;
	*) acx_getkw_dir="-L$with_getkw_libdir" 
	   LDFLAGS="$LDFLAGS $acx_getkw_dir" 
	;;
esac

case $with_getkw in
	yes | "") acx_getkw_libs=getkw;;
	no) acx_getkw_ok=disable ;;
	*) acx_getkw_libs="$with_getkw" ;;
esac

#
# Test if getkw_m.mod is available
#
ACX_F90_MODULE([getkw_m],[acx_getkw_mod_ok=yes],[acx_getkw_mod_ok=no])
AC_LANG_PUSH([Fortran 90])
#
# Test if getkw is linkable
#
# If --with-getkw is defined, then look for THIS AND ONLY THIS lib
if test "x$with_getkw" != x; then
if test $acx_getkw_mod_ok = yes; then
	AC_CHECK_LIB($acx_getkw_libs, end_parse, 
	[acx_getkw_ok=yes], [acx_getkw_ok=specific])
fi; fi

if test $acx_getkw_mod_ok = yes; then
  AC_CHECK_LIB(getkw, end_parse, [acx_getkw_ok=yes;
      acx_getkw_libs=-lgetkw], [acx_getkw_ok=no])
fi
AC_LANG_POP([Fortran 90])

LIBS="$acx_getkw_libs $acx_getkw_save_LIBS"
FFLAGS="$acx_getkw_includes $acx_getkw_save_FFLAGS"
CFLAGS="$acx_getkw_includes $acx_getkw_save_CFLAGS"
LDFLAGS="$acx_getkw_dir $acx_getkw_save_LDFLAGS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test $acx_getkw_ok = no; then
		if test "x$acx_getkw_ok" = "xdisable"; then
			acx_getkw_ok=no
		else
			acx_getkw_ok=no
			$2
		fi
else
		ifelse([$1],,AC_DEFINE(HAVE_LIBGETKW,1,
		[Define if you have the getkw library.]),[$1])
		:
fi
]) dnl ACX_GETKW
