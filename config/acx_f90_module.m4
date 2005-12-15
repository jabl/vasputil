dnl @synopsis ACX_F90_MODULE([MODULE,[ACTION-IF-FOUND, ACTION-IF-NOT-FOUND]])
dnl
dnl Test if a F90 module is available
dnl
dnl @category Fortran90/95
dnl @author Jonas Juselius <jonas@iki.if>
dnl @version 03.05.2003
dnl @license AllPermissive

AC_DEFUN([ACX_F90_MODULE],[
AC_PREREQ(2.57)
AC_REQUIRE([AC_PROG_F90])
m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

AC_LANG_PUSH([Fortran 90])

AC_CACHE_CHECK([for F90 module $1], [acx_cv_f90_mod_$1],[
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [use $1])], 
[acx_cv_f90_mod_$1=yes], [acx_cv_f90_mod_$1=no])])

AC_LANG_POP([Fortran 90])

if test "$acx_cv_f90_mod_$1" = "yes"; then
	$2
	:
else
	$3
	:
fi
]) # ACX_F90_MODULE
