      SUBROUTINE DGESV_F95( A, B, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!     .. "Use Statements" ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO
      USE F77_LAPACK, ONLY: GESV_F77 => LA_GESV
!     .. "Implicit Statement" ..
      IMPLICIT NONE
!     .. "Scalar Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
      REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_GESV computes the solution to a real or complex linear system of
! equations A*X = B, where A is a square matrix and X and B are 
! rectangular matrices or vectors. Gaussian elimination with row 
! interchanges is used to factor A as A = P*L*U , where P is a permutation
! matrix, L is unit lower triangular, and U is upper triangular. The 
! factored form of A is then used to solve the above system.
! 
! =========
! 
!       SUBROUTINE LA_GESV( A, B, IPIV=ipiv, INFO=info )
!          <type>(<wp>), INTENT(INOUT) :: A(:,:), <rhs>
!          INTEGER, INTENT(OUT), OPTIONAL :: IPIV(:)
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
!          <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX square array, shape (:,:).
!        On entry, the matrix A.
!        On exit, the factors L and U from the factorization A = P*L*U; 
!        the unit diagonal elements of L are not stored.
! B      (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(B,1) = size(A,1) or shape (:) with size(B) = size(A,1).
!        On entry, the matrix B.
!        On exit, the solution matrix X .
! IPIV   Optional (output) INTEGER array, shape (:) with size(IPIV) = 
!        size(A,1).
!        The pivot indices that define the permutation matrix P; row i of
!        the matrix was interchanged with row IPIV(i).
! INFO   Optional (output) INTEGER
!        = 0 : successful exit.
!        < 0 : if INFO = -i, the i-th argument has an illegal value.
!        > 0 : if INFO = i, then U(i,i) = 0. The factorization has been 
!        completed, but the factor U is singular, so the solution could
!        not be computed.
!        If INFO is not present and an error occurs, then the program is 
!        terminated with an error message.
!----------------------------------------------------------------------
!     .. "Parameters" ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GESV'
!     .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV
!     .. "Local Pointers" ..
      INTEGER, POINTER :: LPIV(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. "Executable Statements" ..
      LINFO = 0; ISTAT = 0
      IF( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = SIZE(A,1)
      END IF
!     .. "Test the arguments" ..
      IF( SIZE( A, 2 ) /= SIZE(A,1) .OR. SIZE(A,1) < 0 ) THEN
         LINFO = -1
      ELSE IF( SIZE( B, 1 ) /= SIZE(A,1) .OR. SIZE(B,2) < 0 ) THEN
         LINFO = -2
      ELSE IF( SIPIV /= SIZE(A,1) )THEN
            LINFO = -3
      ELSE IF ( SIZE(A,1) > 0 ) THEN
         IF( PRESENT(IPIV) )THEN
            LPIV => IPIV
         ELSE
            ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT )
         END IF
         IF ( ISTAT == 0 ) THEN
!        .. "Call LAPACK77 routine" ..
            CALL GESV_F77( SIZE(A,1), SIZE(B,2), A, MAX(1,SIZE(A,1)), LPIV, B, MAX(1,SIZE(A,1)), &
                           LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT.PRESENT(IPIV) )THEN
            DEALLOCATE(LPIV, STAT = ISTAT1 )
         END IF
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
      END SUBROUTINE DGESV_F95
