AC_DEFUN([ACX_SUBST_FCFLAGS],[

# Fortran flags
case $FC in
gfortran)
	fopt="-O3 -funroll-loops"
	fprof="-pg"
	fdebug="-O0 -g"
	frange="-frange-check"
	fstd="-std=f2003 -Wall -pedantic"
	fldflags=""
	;;
pathf90|pathf95)
	fopt="-Ofast"
	fprof="-pg"
	fdebug="-O0 -g -fno-unsafe-math-optimizations"
	frange="-C"
	fstd="-ansi -fullwarn -Wall"
	fldflags=""
	;;
pgf90|pgf95)
	fopt="-O3 -Bstatic -fastsse -Mipa=fast -Mpreprocess -Minform,warn"
	fprof="-pg"
	fdebug="-O0 -g -Ktrap=fp"
	frange="-C -Mchkfpstk -Mchkptr"
	fstd="-Mstandard"
	fldflags=""
	;;
ifc)
	fopt="-O2 -cm -w95 -fpp"
	fprof="-qp"
	fdebug="-O0 -g"
	frange="-C"
	fldflags="-Vaxlib"
	;;
ifort)
	fopt="-O2 -cm -w95 -fpp"
	fprof="-qp"
	fdebug="-g"
	frange="-C -gen-interface -warn interface -traceback"
	fstd="-warn all -std03"
	fldflags=""
	;;
lf95)
	fopt="--o2 --ntrace --nap --nchk --ng --nsav --prefetch 2 --tpp"
	fprof="-qp"
	fdebug="-g --ap --chkglobal --lst --trace --pca --xref --trap"
	frange="-chkglobal"
	fsd="--f95 --info --wo --warn"
	fldflags="--staticlink"
	;;
xlf95|xlf90)
	fsuffix="-qsuffix=f=f90 -qsuffix=cpp=F90"
        fopt="$fsuffix -O3 -qstrict -qextname"
        fprof="$fsuffix -q"
        fdebug="$fsuffix -g"
        frange="$fsuffix -C"
	fstd="$fsuffix -qlanglvl=95pure"
        fldflags=""
        ;;
f90|fort)
        fopt="-cpp -fast -align dcommons -fpe1"
        fprof="-pg"
        fdebug="-g"
        frange="-C"
        fldflags=""
        ;;
*)
	fopt="-O2"
	finclude="-I/usr/local/include"
	fprof="-q"
	fdebug="-O0 -g"
	frange="-C"
	fldflags=""
	;;
esac

case "${USE_DEBUG_FLAGS}" in
yes)
	fcflags="${fdebug} ${frange} ${fstd}"
	;;
no)
	fcflags="${fopt}"
	;;
*)
	fcflags="${fdebug} ${frange} ${fstd}"
	;;	
esac


AC_SUBST(AM_FCFLAGS, $fcflags)
AC_SUBST(AM_LDFLAGS, $fldflags)

]) #ACX_SUBST_FCFLAGS
