dnl @synopsis ACX_F90_MODULE_NAMING
dnl
dnl Find Fortran 90 modules file extension and case The module extension is
dnl stored in the cached variable acx_cv_f90_modext, or "unknown" if the
dnl extension cannot be found. The file case is stored in acx_cv_f90_modcase
dnl as 'module' or 'MODULE' or 'unknown'. Subsitues F90_MODEXT and F90_MODCASE
dnl on output.
dnl
dnl @category Fortran90/95
dnl @author Jonas Juselius <jonas@iki.if>
dnl @version 03.05.2003
dnl @license AllPermissive

AC_DEFUN([ACX_F90_MODULE_NAMING],[
AC_MSG_CHECKING([Fortran 90 module naming])
AC_CACHE_VAL([acx_cv_f90_modext],
[AC_LANG_PUSH([Fortran 90])
AC_COMPILE_IFELSE([
module conftest
contains
   subroutine test
      return
   end subroutine test
end module conftest
],
  [
for ext in mod MOD d D; do 
  if test -f conftest.$ext; then 
    acx_cv_f90_modext=$ext
    acx_cv_f90_modcase=module
    break
  fi
  if test -f CONFTEST.$ext; then 
    acx_cv_f90_modext=$ext
    acx_cv_f90_modcase=MODULE
  break
  fi
done 
if test x$acx_cv_f90_modext = x ; then
   acx_cv_f90_modext=unknown
   acx_cv_f90_modcase=unknown
else 
   F90_MODEXT=$acx_cv_f90_modext
   F90_MODCASE=$acx_cv_f90_modcase
fi], [acx_cv_f90_modext=unknown])

if test $acx_cv_f90_modcase = module; then
   rm -f conftest.$acx_cv_f90_modext
else
   rm -f CONFTEST.$acx_cv_f90_modext
fi
AC_MSG_RESULT([$acx_cv_f90_modcase.$acx_cv_f90_modext]) 
AC_SUBST(F90_MODEXT)
AC_SUBST(F90_MODCASE)
AC_LANG_POP([Fortran 90])
])])
