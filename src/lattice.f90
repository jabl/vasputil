!****h* vasputil/lattice
! COPYRIGHT
!  Copyright (c) 2004, 2005, 2008 Janne Blomqvist

!  This file is part of Vasputil.

!  Vasputil is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3 of the License, or
!  (at your option) any later version.

!  Vasputil is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with vasputil.  If not, see <http://www.gnu.org/licenses/>.

! PURPOSE
! This module defines a lattice type and procedures for manipulating 
! lattices.
!****

module lattice
  use conf
  use geometry
  use mathconst
!  use la_precision, only: wp => dp
  use f95_lapack

  implicit none

  !****t* lattice/latt
  ! PURPOSE
  ! Type for lattices. Contains the lattice vectors, lattice constant and
  ! volume.
  !****
  type latt
     ! lattice translation vectors for direct (t) and reciprocal (g) lattice.
     ! Vectors are stored in 3x3 matrices as COLUMN vectors
     real(wp), dimension(3,3) :: t, g 
     ! real(wp), dimension(:,:), allocatable :: d ! basis vectors.
     ! Lattice constant and volume of cell
     real(wp) :: a, omega
     ! logical :: cartesian = .FALSE. , reciprocal = .FALSE.     
  end type latt

  contains

    !****f* lattice/latt_init
    ! PURPOSE
    ! Call this subroutine to calculate the reciprocal lattice vectors
    ! and the volume.
    !****
    subroutine latt_init(lattice)
      type(latt), intent(inout) :: lattice
      real(wp) :: rvol
      lattice%omega = dot_product(lattice%t(:,1), lattice%t(:,2) .cross. lattice%t(:,3))
      rvol = 2*pi_Value/lattice%omega
      lattice%g(:,1) = rvol * lattice%t(:,2) .cross. lattice%t(:,3)
      lattice%g(:,2) = rvol * lattice%t(:,3) .cross. lattice%t(:,1)
      lattice%g(:,3) = rvol * lattice%t(:,1) .cross. lattice%t(:,2)

      !print *, 'Reciprocal lattice vectors are: '
      !print *, 'g_1 = ', lattice%g(:,1)
      !print *, 'g_2 = ', lattice%g(:,2)
      !print *, 'g_3 = ', lattice%g(:,3)
    end subroutine latt_init


    !****f* lattice/direct2recip
    ! PURPOSE
    ! Convert a set of points from the direct to the reciprocal lattice.
    ! This would be more efficient by using a Fast Fourier transform, but
    ! this is simple and fast enough.
    !****
    subroutine direct2recip(lattice, coords)
      type(latt), intent(inout) :: lattice
      ! coords should be of size 3xnumAtoms
      real(wp), dimension(:,:), intent(inout) :: coords
      ! We don't want to change the lattice, so we make a copy of it
      real(wp), dimension(3,3) :: rlattvect
      call latt_init(lattice)
      rlattvect = lattice%g

      ! Call LAPACK95 routine to solve Ax=b
      ! print *, 'Calling lapack!'
      call la_gesv (rlattvect, coords)

      ! Multiply by sqrt(3) so that absolute values are comparable
      ! to VASP k-point output list (in OUTCAR).
      coords = coords*sqrt(3.0_wp)
    end subroutine direct2recip


    !****f* lattice/cart2Direct
    ! PURPOSE
    ! Change coordinates from cartesian to direct. This is done by
    ! solving Ax=b where A is the lattice vectors, x are the coordinates
    ! in direct space and b are the coordinates in cartesian space.
    !****
    subroutine cart2Direct(lattice, coords)
      type(latt), intent(in) :: lattice
      real(wp), dimension(:,:), intent(inout) :: coords
      real(wp), dimension(3,3) :: lattvect
      
      lattvect = lattice%t
      
      call la_gesv (lattvect, coords)
    end subroutine cart2Direct


    !****f* lattice/read_lattice
    ! PURPOSE
    ! Read lattice constant and vectors from stdin.
    !****
    subroutine read_lattice (lattice)
      type(latt), intent(out) :: lattice
      integer :: i
      print *, 'Enter the lattice constant:'
      read (*, *) lattice%a
      print *, 'Enter the lattice vectors:'
      do i = 1, 3
         read (*, *) lattice%t(1,i), lattice%t(2,i), &
              lattice%t(3,i)
      end do
      print *, 'read lattice stuff'
      if (is_lin_dependent (lattice%t)) then
         call error_stop ("ERROR: Lattice vectors are not linearly independent!")
      end if
      call latt_init (lattice)
    end subroutine read_lattice

end module lattice
