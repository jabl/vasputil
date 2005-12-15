dnl @synopsis ACX_F90_WITH_LIBDIR([])
dnl
dnl Simply set FFLAGS and FLDFLAGS
dnl
dnl @category Fortran90/95
dnl @author Jonas Juselius <jonas@iki.if>
dnl @version 03.05.2003
dnl @license AllPermissive

AC_DEFUN([ACX_F90_WITH_LIBDIR],[
AC_PREREQ(2.57)

AC_ARG_WITH([f90_libdir], AC_HELP_STRING([--with-f90-libdir=DIR],
	[Base path to installed f90 libraries)]))

AC_ARG_WITH([f90_moddir], 
  AC_HELP_STRING([--with-f90-modules=DIR],
  [Directory where F90 modules and include files can be found]))
  
AC_ARG_WITH([f90_libs], 
  AC_HELP_STRING([--with-f90-libs=DIR],
  [Directory where F90 library files can be found]))
  
case $with_f90_libdir in
	  yes | no | "") ;;
	  *) LDFLAGS="$LDFLAGS -L$with_f90_libdir/lib" 
		 FFLAGS="$FFLAGS -I$with_f90_libdir/include"
		 CFLAGS="$CFLAGS -I$with_f90_libdir/include"
	  ;;
esac

case $with_f90_moddir in
	yes | no | "") ;;
	*) FFLAGS="$FFLAGS -I$with_f90_moddir"
	   CFLAGS="$CFLAGS -I$with_f90_moddir"
	 ;;
esac

case $with_f90_libs in
	yes | no | "") ;;
	*) LDFLAGS="$LDFLAGS -L$with_f90_libs" 
	;;
esac
]) # ACX_F90_MODULE
