!****h* vasputil/conf
! PURPOSE
! Global configuration data for the application
!****
module conf

  ! Working precision is IEEE double precision
  use kind_params, only: wp => dp, i8b

  implicit none

  ! Change this to compile with(out) debug prints.
  logical, parameter :: debug = .false.

  contains


    !****f* conf/error_stop
    ! PURPOSE
    !   Print an error message and stop the program.
    !****
    subroutine error_stop(message)
      character(len=*), intent(in), optional :: message
      if (present(message)) then
         print *, message
      end if
      print *, 'Error encountered, stopping.'
      stop 1
    end subroutine error_stop


    !****f* conf/error_msg
    ! PURPOSE
    !  Depending on the presence of the optional parameter,
    !  set the status flag or call error_stop
    !****
    subroutine error_msg (flagval, message, status)
      character(len=*), intent(in), optional :: message
      integer, intent(in) :: flagval
      integer, intent(out), optional :: status

      if (present (status)) then
         status = flagval
      else
         call error_stop (message)
      end if
    end subroutine error_msg


end module conf
