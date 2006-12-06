dnl @synopsis ACX_FC_SEARCH_PATH([])
dnl
dnl Simply set FCFLAGS and LDFLAGS
dnl The idea here is that Fortran development libraries, together with headers
dnl and module files, should be installed in a standard way so that they are
dnl easy to locate.
dnl 
dnl @category Fortran90/95
dnl @author Jonas Juselius <jonas@iki.if>
dnl @version 13.05.2006
dnl @license AllPermissive

AC_DEFUN([ACX_FC_SEARCH_PATH],[
AC_PREREQ(2.59)

AC_ARG_WITH([fc_search_path], AC_HELP_STRING([--with-fc-search-path=DIR],
	[Base path to installed Fortan libraries)]))

dnl AC_ARG_WITH([fmoddir], 
dnl   AC_HELP_STRING([--with-fc-modules-path=DIR],
dnl   [Directory where Fortran module and include files can be found]))
dnl   
dnl AC_ARG_WITH([flibs], 
dnl   AC_HELP_STRING([--with-fc-lib-path=DIR],
dnl   [Directory where Fortran library files can be found]))

case $with_fc_search_path in
	  yes | no | "") ;;
	  *) LDFLAGS="$LDFLAGS -L$with_fc_search_path/lib" 
		 FCFLAGS="$FCFLAGS -I$with_fc_search_path/include"
	  ;;
esac

dnl case $with_fmoddir in
dnl     yes | no | "") ;;
dnl     *) FCFLAGS="$FCFLAGS -I$with_fmoddir"
dnl      ;;
dnl esac

dnl case $with_flibs in
dnl     yes | no | "") ;;
dnl     *) LDFLAGS="$LDFLAGS -L$with_flibs" 
dnl     ;;
dnl esac
]) # ACX_FC_SEARCH_PATH
