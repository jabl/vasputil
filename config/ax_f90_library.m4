dnl @synopsis AX_F90_LIBRARY(LIBRARY, LIB-REGEXP, FUNCTION-BODY [, SEARCH-PATH [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Set up the compiler flags to link a given fortran 90 library
dnl LIBRARY is the name of the library. LIB-REGEXP is a regular
dnl expression (used by find) matched by the filename of the library,
dnl this is useful either if the library filename does not follow the
dnl traditional libxxx.a or libxxx.so pattern, or if some specific
dnl information is embedded into the name, like compiler used,
dnl debugging status ...). FUNCTION-BODY is the body of a function
dnl (including the 'use' statements and the call to a function defined
dnl by the library) SEARCH-PATH is a colon-separated list of
dnl directories that will be used as the base directoris for 'find' to
dnl look for the library file. If empty, the search path will be
dnl composed of $prefix/lib, $ac_default_prefix/lib, and
dnl $LD_LIBRARY_PATH. Two output variables named F90_LDFLAGS_xxx and
dnl F90_LIBS_xxx will be set up with the proper flag for substitution
dnl in Makefiles (xxx is built from the first argument, with autoconf
dnl traditional escapes).
dnl
dnl @category Fortran
dnl @author Luc Maisonobe <luc@spaceroots.org>
dnl @version 2005-01-14
dnl @license AllPermissive

AC_DEFUN([AX_F90_LIBRARY],[
AS_VAR_PUSHDEF([ax_ldflags],[ax_f90_ldflags_$1])
AS_VAR_PUSHDEF([ax_libs],[ax_f90_libs_$1])
AC_MSG_CHECKING([$1 fortran 90 library])
AC_LANG_PUSH(Fortran)
AS_VAR_SET([ax_ldflags],"")
AS_VAR_SET([ax_f90_library_$1_ok],"no")
AS_VAR_SET([ax_libs],"not found")

if test "x$4" = x ; then
  ax_search=""
else 
  ax_search="$4"
fi

if test "$prefix" = "NONE"; then
  ax_search="$ax_search:$ac_default_prefix/lib"
else
  ax_search="$ax_search:$prefix/lib:$ac_default_prefix/lib"
fi

for ax_base in "" `echo $LD_LIBRARY_PATH | tr ':' '\012'` ; do
  if test "x$ax_base" != x ; then
    changequote(,)dnl
    ax_base=`echo $ax_base | sed 's,/[^/]*$,,'`
    changequote([,])dnl
    ax_search="${ax_search}:${ax_base}"
  fi
done

ax_save_LDFLAGS="$LDFLAGS"
ax_save_LIBS="$LIBS"

for ax_base in `echo $ax_search | tr ':' '\012'` ; do
 if test "AS_VAR_GET(ax_libs)" = "not found" ; then
   for ax_lib in "" `find $ax_base -follow -name '$2' -print` ; do
     if test "x$ax_lib" != x ; then
       changequote(,)dnl
       ax_dir=`echo $ax_lib | sed 's,/[^/]*$,,'`
       ax_lib=`echo $ax_lib | sed 's,.*/\([^/]*\)$,\1,'`
       changequote([,])dnl
       case "$ax_lib" in
         lib*)
           changequote(,)dnl
           ax_lib="`echo $ax_lib | sed 's,lib\(.*\)\.[^.]*$,\1,'`"
           changequote([,])dnl
           AS_VAR_SET([ax_ldflags],"-L$ax_dir")
           AS_VAR_SET([ax_libs],"-l$ax_lib")
           ;;
         *)
           AS_VAR_SET([ax_ldflags],"")
           AS_VAR_SET(ax_libs,"$ax_lib")
           ;;
       esac
       LDFLAGS="$ax_save_LDFLAGS AS_VAR_GET(ax_ldflags)"
       LIBS="AS_VAR_GET(ax_libs) $ax_save_LIBS"
       AC_LINK_IFELSE([program conftest_program
$3
          end program conftest_program
         ],[],[AS_VAR_SET(ax_ldflags,"")
          AS_VAR_SET(ax_libs,"not found")
         ])
     fi
   done
 fi
done
AC_LANG_POP(Fortran)
AC_MSG_RESULT([AS_VAR_GET(ax_ldflags) AS_VAR_GET(ax_libs)])
if test "AS_VAR_GET(ax_libs)" = "not found"; then
 AS_TR_SH(F90_LDFLAGS_$1)=""
 AS_TR_SH(F90_LIBS_$1)=""
 $5
else
 AS_TR_SH(F90_LDFLAGS_$1)=AS_VAR_GET(ax_ldflags)
 AS_TR_SH(F90_LIBS_$1)=AS_VAR_GET(ax_libs)
 AS_VAR_SET([ax_f90_library_$1_ok],"yes")
 $6
fi
AC_SUBST(AS_TR_SH(F90_LDFLAGS_$1))
AC_SUBST(AS_TR_SH(F90_LIBS_$1))
AS_VAR_POPDEF([ax_libs])
AS_VAR_POPDEF([ax_ldflags])
])
