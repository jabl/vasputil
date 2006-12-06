dnl @synopsis AX_F90_MODULE(MODULE, MODULE-REGEXP, FUNCTION-BODY [, SEARCH-PATH [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Set up the compiler flags to use a given fortran 90 module MODULE
dnl is the name of the module. MODULE-REGEXP is a regular expression
dnl (used by find) matched by the filename of the module. FUNCTION-BODY
dnl is the body of a function (including the 'use' statement and the
dnl call to a function defined by the module) SEARCH-PATH is a
dnl colon-separated list of directories that will be recursively
dnl searched for modules files. If empty, the search path will be
dnl composed of $prefix, $ac_default_prefix, and all directories
dnl exactly one level *above* the directories in $LD_LIBRARY_PATH (the
dnl rationale is that when libraries are put in /some/path/lib, the
dnl modules are often put in a directory like /some/path/include or
dnl /some/path/mod or something similar). An output variable named
dnl F90_MODULE_xxx will be set up with the proper flag for substitution
dnl in Makefiles (xxx is built from the first argument, with autoconf
dnl traditional escapes).
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-01-14
dnl @license AllPermissive

AC_DEFUN([AX_F90_MODULE],[
 AC_REQUIRE([AX_F90_MODULE_FLAG])
 AX_F90_INTERNAL_HEADMOD([$1 fortran 90 module],[$2],"$ax_f90_modflag",
                         [$3],AS_TR_SH(F90_MODULE_$1),[$4],[$5],[$6])
])
