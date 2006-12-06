dnl @synopsis AX_F90_LIBRARY_SETUP(LIBRARY, HEADER-REGEXP, MODULE-REGEXP, LIB-REGEXP, FUNCTION-BODY [, PATH])
dnl
dnl Convenience macro to set up a fortran 90 library in a simplified
dnl way. LIBRARY is the name of the library. HEADER-REGEXP is a regular
dnl expression (used by find) matched by the header file to look for
dnl (may be empty). MODULE-REGEXP is a regular expression (used by
dnl find) matched by the filename of the module (may be empty).
dnl LIB-REGEXP is a regular expression (used by find) matched by the
dnl filename of the library, this is useful either if the library
dnl filename does not follow the traditional libxxx.a or libxxx.so
dnl pattern, or if some specific information is embedded into the name,
dnl like compiler used, debugging status ...). FUNCTION-BODY is the
dnl body of a function (including the 'use' statements and the call to
dnl a function defined by the library).
dnl
dnl This macro is a simple wrapper around AX_F90_MODULE and
dnl AX_F90_LIBRARY that uses the parameters provided by the end user
dnl through a --with-xxx option to set up the search path. Both a
dnl module and a library will be tested, the same path will be used for
dnl both tests, so the path must be set up with a common parent
dnl directory of both the library file and the module file. The macro
dnl also automatically updates the FCFLAGS, LDFLAGS and LIBS variables
dnl in addition to providing the F90_HEADER_xxx, F90_MODULE_xxx,
dnl F90_LDFLAGS_xxx and F90_LIBS_xxx output variables.
dnl
dnl Example: suppose you have
dnl /home/nostradamus/esoteric/lib/libalchemy.a and
dnl /home/nostradamus/esoteric/mod/alchemy.mod which provides a
dnl function transmute_into_gold, you can use the following in you
dnl configure.ac:
dnl
dnl   AX_F90_MODULE_EXTENSION
dnl   if test x$ax_f90_modext = xunknown ; then
dnl     AC_MSG_ERROR([unable to find f90 modules extension])
dnl   fi
dnl   AX_F90_LIBRARY_SETUP(alchemy,[],alchemy.$ax_f90_modext,libalchemy*,[
dnl      use alchemy
dnl      call transmute_into_gold('lead')
dnl     ])
dnl
dnl and the user could configure your package using a command like
dnl this:
dnl
dnl   ./configure --with-alchemy=$HOME/esoteric
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-01-14
dnl @license AllPermissive

AC_DEFUN([AX_F90_LIBRARY_SETUP],[
 AC_ARG_WITH([$1],AC_HELP_STRING([--with-$1=PATH], 
 [specify search path form $1 module and library]),
   [if test x${withval} = xno ; then
     AC_MSG_WARN([$1 disabled at user option])
    fi],[withval="$6"])
 if test x$2 != x ; then
   AX_F90_HEADER([$1],[$2],[$5],$withval,[
      FCFLAGS="$FCFLAGS $AS_TR_SH(F90_HEADER_$1)"
     ],[])
 fi
 if test x$3 != x ; then
   AX_F90_MODULE([$1],[$3],[$5],$withval,[
      FCFLAGS="$FCFLAGS $AS_TR_SH(F90_MODULE_$1)"
     ],[])
 fi
 AX_F90_LIBRARY([$1],[$4],[$5],$withval,[
    LDFLAGS="$LDFLAGS $AS_TR_SH(F90_LDFLAGS_$1)"
    LIBS="$AS_TR_SH(F90_LIBS_$1) $LIBS"
   ],[])
])
