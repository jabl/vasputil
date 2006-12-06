dnl
dnl @synopsis ACX_F90_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Written by Jonas Juselius <jonas@iki.fi>
dnl

AC_DEFUN([ACX_F90_MPI], [
AC_PREREQ(2.59)

ACX_FC_LIBRARY_SETUP([mpi_f90],[mpi],[libmpi_f90*],[
use mpi_kinds 
use mpi 
integer :: i 
call mpi_init(i) 
],[],[$1],[$2])

]) dnl ACX_F90_MPI
