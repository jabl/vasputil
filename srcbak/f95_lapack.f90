!  Copyright (c) 2004 Janne Blomqvist

! This module blatantly steals some stuff from lapack95 (Fortran 95
! interface to lapack). If you have the real thing, you can replace
! this file with lapack95.

module f95_lapack
  implicit none

  INTERFACE LA_GESV

     SUBROUTINE DGESV_F95( A, B, IPIV, INFO )
       USE LA_PRECISION, ONLY: WP => DP
       INTEGER, INTENT(OUT), OPTIONAL :: INFO
       INTEGER, INTENT(OUT), OPTIONAL :: IPIV(:)
       REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
     END SUBROUTINE DGESV_F95

  END INTERFACE

  INTERFACE LA_GESDD

     SUBROUTINE DGESDD_F95(A, S, U, VT, WW, JOB, INFO )
       USE LA_PRECISION, ONLY: WP => DP
       CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: JOB
       INTEGER, INTENT(OUT), OPTIONAL :: INFO
       REAL(WP), INTENT(INOUT) :: A(:,:)
       REAL(WP), INTENT(OUT) :: S(:)
       REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: WW(:)
       REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: U(:,:), VT(:,:)
     END SUBROUTINE DGESDD_F95

  END INTERFACE

end module f95_lapack
