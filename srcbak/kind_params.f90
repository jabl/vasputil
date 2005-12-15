!****h* /kind_params
! PURPOSE
! Kind type parameters.
!****
module kind_params

  implicit none

  ! 1, 2, 4 and 8 byte integers.
  integer, parameter :: i1b = selected_int_kind(2)
  integer, parameter :: i2b = selected_int_kind(4)
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i8b = selected_int_kind(18)

  ! IEEE Single precision floating point.
  integer, parameter :: sp = selected_real_kind(6)

  !  integer, parameter :: dp = kind(1.0d0) ! Double precision
  
  ! At least IEEE 754 double precision.
  integer, parameter :: dp = selected_real_kind(15, 307) 
  
  ! quad precision, fall back on double if quad is not available.
  integer, parameter :: qp_preferred = selected_real_kind(30,1000)
  ! Portland compiler doesn't understand the following although it is
  ! valid Fortran 95.
  !integer, parameter :: qp = (1+sign(1,qp_preferred))/2*qp_preferred+ &
  !                            (1-sign(1,qp_preferred))/2*dp

end module kind_params
