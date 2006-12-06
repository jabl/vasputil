dnl
dnl @synopsis ACX_FC_LIBRARY_SETUP()
dnl
dnl Written by Jonas Juselius <jonas@iki.fi>
dnl

AC_DEFUN([ACX_FC_LIBRARY_SETUP], [
AC_PREREQ(2.59)

if test "x$ax_f90_modext" = "xunknown" ; then
     AC_MSG_ERROR([unable to find f90 modules extension])
fi

AS_VAR_SET([acx_$1_ok],["no"])

AC_ARG_WITH([$1], AC_HELP_STRING([--with-$1=DIR],
[specify search path form $1 module and library]),
[case $withval in
   yes|"") 
       withval=$FC_SEARCH_PATH ;;
   no)
       AC_MSG_WARN([blas95 disabled at user option])
	   withval=""
       acx_blas95_ok="disabled"
   ;;
    esac],
[withval="$FC_SEARCH_PATH"])

if test x"AS_VAR_GET(acx_$1_ok)" != xdisabled; then

ax_blas95_save_LIBS="$LIBS"
AS_VAR_SET([acx_$1_save_LIBS],["$LIBS"])
LIBS="$LIBS $5"

AX_F90_MODULE([$1],[$2.$ax_f90_modext],[$4],$withval,
[FCFLAGS="$FCFLAGS $AS_TR_SH(F90_MODULE_$1)"],
[AS_VAR_SET([acx_$1_ok],["yes"])])


AX_F90_LIBRARY([$1],[$3],[$4],$withval,
[LDFLAGS="$LDFLAGS $AS_TR_SH(F90_LDFLAGS_$1)" 
LIBS="$AS_TR_SH(F90_LIBS_$1) $LIBS"],
[AS_VAR_SET([acx_$1_ok],["yes"])])

#LIBS="$ax_blas95_save_LIBS $F90_LIBS_blas95"
LIBS="AS_VAR_GET(acx_$1_save_LIBS)"
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"AS_VAR_GET(acx_$1_ok)" = xyes; then
		ifelse([$1],,i
		AC_DEFINE(AS_TR_SH(HAVE_$1),1,[Define if you have $1 library.]),[$6])
		:
else
        if test x"AS_VAR_GET(acx_$1_ok)" = xdisabled; then
          AS_VAR_SET([acx_$1_ok],["no"])
		  :
		else
          AS_VAR_SET([acx_$1_ok],["no"])
		  $7
		fi
fi

]) dnl ACX_FC_LIBRARY_SETUP
