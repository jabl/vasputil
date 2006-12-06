dnl @synopsis AX_F90_HEADER(HEADER, HEADER-REGEXP, FUNCTION-BODY [, SEARCH-PATH [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Set up the compiler flags to use a given fortran 90 header. HEADER
dnl is the name of the header. HEADER-REGEXP is a regular expression
dnl (used by find) matched by the filename of the header. FUNCTION-BODY
dnl is the body of a function (including the 'use' statement and the
dnl call to a function defined by the module) SEARCH-PATH is a
dnl colon-separated list of directories that will be recursively
dnl searched for header files. If empty, the search path will be
dnl composed of $prefix, $ac_default_prefix, and all directories
dnl exactly one level *above* the directories in $LD_LIBRARY_PATH (the
dnl rationale is that when libraries are put in /some/path/lib, the
dnl headers are often put in a directory like /some/path/include). An
dnl output variable named F90_HEADER_xxx will be set up with the proper
dnl flag for substitution in Makefiles (xxx is built from the first
dnl argument, with autoconf traditional escapes).
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-01-14
dnl @license AllPermissive

AC_DEFUN([AX_F90_HEADER],[
 AX_F90_INTERNAL_HEADMOD([$1 fortran 90 header],[$2],-I,
                         [$3],AS_TR_SH(F90_HEADER_$1),[$4],[$5],[$6])
])
