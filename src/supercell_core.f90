!****h* vasputil/supercell_core
! NAME
!   supercell_core
! COPYRIGHT
!  Copyright (c) 2004, 2005, 2006, 2008 Janne Blomqvist

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
! This module defines types for supercells and procedures for manipulating 
! these.
!****
module supercell_core
  use conf
  use geometry
  use lattice

  implicit none


  !****t* supercell_core/atom
  ! PURPOSE
  !   Type for an atom, containing the atomic sumbol and the 3D selective flag.
  !****
  type atom
     character(len = 2) :: symbol = " H"
     logical, dimension(3) :: selective = .true.
  end type atom


  !****t* supercell_core/supercell
  ! PURPOSE
  !   A supercell type. Contains an array of atom types and their coordinates.
  !****
  type supercell
     type(latt) :: lattice
     ! Fortran 95 doesn't allow allocatable components
     ! in derived types. Thus we use pointers.
     real(wp), dimension(:, :), pointer :: atomCoords => NULL()
     type(atom), dimension(:), pointer :: atoms => NULL()
     logical :: cartesian = .FALSE. , relative = .TRUE., selective = .false.
  end type supercell

contains


  !****f* supercell_core/direct2Cartesian
  ! PURPOSE
  ! Change supercell coordinates from direct to cartesian.
  !****
  subroutine direct2Cartesian(cell)
    type(supercell), intent(inout) :: cell
    !      print *, 'latt. vect size: ', size(cell%latticeVectors, 1)
    !      print *, 'atomCoord matrix size: ', size(cell%atomCoords, 1)
    if (cell%cartesian) then
       return
    end if
    cell%atomCoords = matmul( cell%lattice%t, &
         cell%atomCoords)
    cell%cartesian = .TRUE.
  end subroutine direct2Cartesian


  !****f* supercell_core/cartesian2Direct
  ! PURPOSE
  ! Change supercell coordinates from cartesian to direct.
  !****
  subroutine cartesian2Direct(cell)
    type(supercell), intent(inout) :: cell
    if (.not. cell%cartesian) then
       return
    end if
    !call error_stop('cartesian2Direct not implemented yet!')
    call cart2Direct(cell%lattice, cell%atomCoords)
    cell%cartesian = .false.
  end subroutine cartesian2Direct


  !****f* supercell_core/rel2Act
  ! PURPOSE
  ! Change supercell coordinates from relative to actual.
  !****
  subroutine rel2Act(cell)
    type(supercell), intent(inout) :: cell
    if (cell%relative) then
       cell%atomCoords = cell%atomCoords * cell%lattice%a
       cell%relative = .FALSE.
    end if
  end subroutine rel2Act


  !****f* supercell_core/act2Rel
  ! PURPOSE
  ! Change supercell coordinates from actual to relative.
  !****
  subroutine act2Rel(cell)
    type(supercell), intent(inout) :: cell
    if (.not. cell%relative) then
       cell%atomCoords = cell%atomCoords / cell%lattice%a
       cell%relative = .TRUE.
    end if
  end subroutine act2Rel


  !****f* supercell_core/normalize
  ! PURPOSE
  ! Normalize the atom coordinates in a cell, i.e. make sure that
  ! the relative coordinates of all atoms are between 0 and 1.
  !****
  subroutine normalize(cell)
    type(supercell), intent(inout) :: cell
    if (cell%cartesian) then
       call cartesian2Direct(cell)
    end if
    if (.not. cell%relative) then
       call act2Rel(cell)
    end if
    cell%atomCoords = mod(cell%atomCoords, 1.0_wp)
    cell%atomCoords = cell%atomCoords + 1.0_wp
    cell%atomCoords = mod(cell%atomCoords, 1.0_wp)
  end subroutine normalize


end module supercell_core
