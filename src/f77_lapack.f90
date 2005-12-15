!  Copyright (c) 2004 Janne Blomqvist

! This module blatantly steals some stuff from lapack95 (Fortran 95
! interface to lapack). If you have the real thing, you can replace
! this file with lapack95.

module f77_lapack

  implicit none

  INTERFACE LA_GESV

     SUBROUTINE DGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
       USE LA_PRECISION, ONLY: WP => DP
       INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
       INTEGER, INTENT(OUT) :: INFO
       INTEGER, INTENT(OUT) :: PIV(*)
       REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
     END SUBROUTINE DGESV

  END INTERFACE


  INTERFACE LA_GESDD


     SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,    &
          &                      WORK, LWORK, IWORK, INFO )
       USE LA_PRECISION, ONLY: WP => DP
       CHARACTER(LEN=1), INTENT(IN) :: JOBZ
       INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
       INTEGER, INTENT(OUT) :: INFO
       REAL(WP), INTENT(OUT) :: S(*)
       REAL(WP), INTENT(INOUT) :: A(LDA,*)
       REAL(WP), INTENT(OUT) :: U(LDU,*), VT(LDVT,*), WORK(*)
       INTEGER :: IWORK(*)
     END SUBROUTINE DGESDD


  END INTERFACE


end module f77_lapack
