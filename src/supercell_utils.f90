!****h* vasputil/supercell_utils
! NAME
!   supercell_utils
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
! This module defines some auxiliary helper procedures for supercells.
!****
module supercell_utils
  use conf
  use geometry
  use lattice
  use supercell_core

  implicit none

contains


  !****f* supercell_core/init_cell
  ! PURPOSE
  ! Constructor for cell types.
  !****
  subroutine init_cell (cell, size, stat)
    type(supercell), intent(inout) :: cell
    integer, intent(in) :: size
    integer, optional, intent(out) :: stat
    ! Clean out old goop, if present
    call destroy_cell (cell, stat)
    if (present (stat)) then
       allocate (cell%atoms(size), cell%atomCoords(3, size), stat=stat)
    else
       allocate (cell%atoms(size), cell%atomCoords(3, size))
    end if
  end subroutine init_cell


  !****f* supercell_core/destroy_cell
  ! PURPOSE
  ! Destructor for cell types.
  !****
  subroutine destroy_cell (cell, stat)
    type(supercell), intent(inout) :: cell
    integer, optional, intent(out) :: stat
    if (associated (cell%atoms) .and. associated (cell%atomCoords)) then
       if (present (stat)) then
          deallocate (cell%atoms, cell%atomCoords, stat=stat)
       else
          deallocate (cell%atoms, cell%atomCoords)
       end if
    else ! This shouldn't happen..
       if (present (stat)) then
          stat = -9999
       end if
    end if
  end subroutine destroy_cell


  !****f* supercell_core/sort_cell
  ! PURPOSE
  ! Sort a cell so that atoms of the same species are next to each other.
  !****
  subroutine sort_cell (cell)
    type(supercell), intent(inout) :: cell

    if (debug) then
       print *, 'qsort: lbound: ', lbound (cell%atoms, dim=1)
       print *, 'qsort: ubound: ', ubound (cell%atoms, dim=1)
    end if

    call qsort (lbound (cell%atoms, dim=1), ubound (cell%atoms, dim=1))

  contains

    !****f* sort_cell/qsort
    ! PURPOSE
    ! The quicksort algorithm, with the middle array
    ! element chosen as pivot.
    !****
    recursive subroutine qsort (begin, end)
      integer, intent(in) :: begin, end
      integer :: ipivot, left, right
      character(len=2) :: apivot

      !        print *, 'qsort called with begin: ', begin, ' end: ', end
      if (end > begin + 1) then
         left = begin
         right = end
         ipivot = (begin + end) / 2
         apivot = adjustl (cell%atoms(ipivot)%symbol)
         !           call swap (begin, ipivot)
         qout: do
            do while (llt (adjustl (cell%atoms(left)%symbol), apivot))
               left = left + 1
               !                 print *, 'increasing left'
            end do
            do while (llt ( apivot, adjustl (cell%atoms(right)%symbol)))
               right = right - 1
               !                 print *, 'decreasing right'
            end do
            if (left < right) then
               call swap (left, right)
               if (ipivot == left) then
                  ipivot = right
               else if (ipivot == right) then
                  ipivot = left
               end if
               left = left + 1
               right = right - 1
            else if (left == right) then
               left = left + 1
               right = right - 1
               exit qout
            else
               exit qout
            end if
         end do qout
         !           print *, 'left: ', left, ' right: ', right
         !           left = left - 1
         !           call swap (begin, left)
         call qsort (begin, right)
         call qsort (left, end)
      else if (end - begin == 1) then
         if (llt (adjustl (cell%atoms(end)%symbol), &
              adjustl (cell%atoms(begin)%symbol))) then
            call swap (begin, end)
         end if
      end if
    end subroutine qsort

    !****f* sort_cell/swap
    ! PURPOSE
    ! Swap two atoms.
    !****
    subroutine swap (i, j)
      integer, intent(in) :: i, j
      type(atom) :: tmpatom
      real(wp) :: tmpcoords(3)

      tmpatom = cell%atoms(i)
      !        if (debug) then
      !           print *, 'tmpatom: ', tmpatom
      !        end if
      tmpcoords = cell%atomCoords(:, i)
      cell%atoms(i) = cell%atoms(j)
      cell%atomCoords(:, i) = cell%atomCoords(:, j)
      cell%atoms(j) = tmpatom
      cell%atomCoords(:, j) = tmpcoords
    end subroutine swap

  end subroutine sort_cell


  !****f* supercell_utils/checkCells
  ! PURPOSE
  ! Check that two supercells are consistent so that e.g. interpolation
  ! is sensible.
  !****
  function checkCells(cell1, cell2, errmsg)
    type(supercell), intent(in) :: cell1, cell2
    character(len=*), intent(out), optional :: errmsg !Error message
    logical :: checkCells

    if (present (errmsg)) then
       errmsg = ''
    end if
    ! We start by assuming they are compatible
    checkCells = .TRUE.
    ! Check that latticeConstants match
    if ( abs(cell1%lattice%a - cell2%lattice%a) > &
         sqrt(epsilon(cell1%lattice%a))) then
       if (present(errmsg)) then
          errmsg = "Lattice constants of cells don't match."
       end if
       checkCells = .FALSE.
       return
    end if
    ! Check lattice vectors.
    if (any( abs(cell1%lattice%t - cell2%lattice%t) > &
         sqrt( epsilon(cell1%lattice%t)) )) then
       if (present(errmsg)) then
          errmsg = "Lattice vectors of cells don't match."
       end if
       checkCells = .FALSE.
       return
    end if
    ! Check number of atoms
    if (size(cell1%atoms) /= size(cell2%atoms)) then
       if (present(errmsg)) then
          errmsg = 'Number of atoms in cells not consistent.'
       end if
       checkCells = .FALSE.
       return
    end if
  end function checkCells


  !****f* supercell_utils/check_nndist
  ! PURPOSE
  ! Check that all atoms in a supercell are at least a specified distance
  ! from each other.
  !****
  subroutine check_nndist(cell, tol)
    type(supercell), intent(inout) :: cell
    real(wp), intent(in), optional :: tol
    real(wp) :: t, dist
    integer :: i,j

    if (present(tol)) then
       t = tol
    else
       t = 1.0_wp ! Minimum 1 Å
    end if

    call direct2Cartesian(cell)
    call rel2Act(cell)

    do j = 1, size(cell%atoms)
       do i = 1, j-1
          dist = norm(cell%atomcoords(:, i) - cell%atomcoords(:, j))
          if (dist < t) then
             print *, 'WARNING: Distance between atoms ', i, ' and ', &
                  j, ' is only ', dist, ' Å!'
          end if
       end do
    end do
  end subroutine check_nndist


end module supercell_utils
