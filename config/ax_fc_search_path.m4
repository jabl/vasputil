dnl @synopsis AX_FC_SEARCH_PATH([])
dnl
dnl Simply set FCFLAGS and LDFLAGS. Additionally set FC_SEARCH_PATH for
dnl nice interplay with the AX_F90_* macros.
dnl The idea here is that Fortran development libraries, together with headers
dnl and module files, should be installed in a standard way so that they are
dnl easy to locate.
dnl 
dnl @category Fortran90/95
dnl @author Jonas Juselius <jonas@iki.if>
dnl @version 13.05.2006
dnl @license AllPermissive

AC_DEFUN([AX_FC_SEARCH_PATH],[
AC_PREREQ(2.59)

AC_ARG_WITH([fc_path], AC_HELP_STRING([--with-fc-path=PATH],
	[Base path to installed Fortran libraries]))

case $with_fc_path in
	  yes | no | "") ;;
	  *) LDFLAGS="$LDFLAGS -L$with_fc_path/lib" 
		 FCFLAGS="$FCFLAGS -I$with_fc_path/include"
		 CFLAGS="$CFLAGS -I$with_fc_path/include"
		 if test "x$FC_SEARCH_PATH" = x; then
		    FC_SEARCH_PATH="$with_fc_path"
		 else
		    FC_SEARCH_PATH="$FC_SEARCH_PATH:$with_fc_path"
		 fi
	  ;;
esac

]) # AX_FC_SEARCH_PATH
