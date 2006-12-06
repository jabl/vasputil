dnl
dnl @synopsis ACX_BLAS95([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Written by Jonas Juselius <jonas@iki.fi>
dnl

AC_DEFUN([ACX_BLAS95], [
AC_PREREQ(2.59)

ACX_FC_LIBRARY_SETUP([blas95],[blas95],[libblas95*],[
use blas95
real(8), dimension(2,2) :: a=0.0, b=0.0, c
call gemm(a,b,c)
],[$BLAS_LIBS],[$1],[$2])
]) dnl ACX_BLAS95
