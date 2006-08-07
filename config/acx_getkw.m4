dnl
dnl @synopsis ACX_GETKW([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl Written by Jonas Juselius <jonas@iki.fi>
dnl

AC_DEFUN([ACX_GETKW], [
AC_PREREQ(2.59)

ACX_FC_LIBRARY_SETUP([getkw],[getkw_m],[libgetkw*],[
use getkw_m
call set_parser_verbose(0)
],[],[$1],[$2])

]) dnl ACX_GETKW
