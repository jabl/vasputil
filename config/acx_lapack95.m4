dnl
dnl @synopsis ACX_LAPACK95([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Written by Jonas Juselius <jonas@iki.fi>
dnl

AC_DEFUN([ACX_LAPACK95], [
AC_PREREQ(2.59)

ACX_FC_LIBRARY_SETUP([lapack95],[lapack95],[liblapack95*],[
use lapack95
real(8), dimension(2,2) :: a=0.0
real(8), dimension(2) :: b=0.0
call syev(a,b)
],[$LAPACK_LIBS $BLAS_LIBS],[$1],[$2])
]) dnl ACX_LAPACK95
