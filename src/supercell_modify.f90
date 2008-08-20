!****h* vasputil/supercell_modify
! NAME
!   supercell_modify
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
!   This module defines procedures that modify a supercell.
!****
module supercell_modify
  use conf
  use supercell_core
  use poscar_io
  use supercell_utils

  implicit none

contains

  !****f* supercell_modify/normalize_POSCAR
  ! PURPOSE
  !   Normalize the coordinates, i.e. make sure that they are between 0 and
  !   1 in direct coordinates, in a POSCAR file.
  !****
  subroutine normalize_POSCAR(infile, outfile)
    character(len=*), intent(in) :: infile, outfile
    type(supercell) :: cell
    call read_POSCAR(cell, infile)
    call normalize(cell)
    call write_POSCAR(cell, outfile)
  end subroutine normalize_POSCAR


  !****f* supercell_modify/unnormalize_POSCAR
  ! PURPOSE
  !   "Unnormalize" a POSCAR file, i.e. try to make the atoms fit in a box.
  ! NOTES
  !   Somewhat ad hoc and semimanual, don't expect magic results.
  !****
  subroutine unnormalize_POSCAR(infile, dir, atoms, outfile)
    character(len=*), intent(in) :: infile, outfile
    integer, intent(in) :: dir, atoms(:)
    type(supercell) :: cell
    call read_POSCAR (cell, infile)
    if (debug) then
       print *, 'Unnormalizing dir: ', dir, ' atoms: ', atoms
    end if
    call unnormalize (cell, dir, atoms)
    call write_POSCAR (cell, outfile)
  end subroutine unnormalize_POSCAR


  !****f* supercell_core/unnormalize
  ! PURPOSE
  ! "Unnormalize" atom coordinates in a cell. Useful when one wants
  ! to make a supercell box-like e.g. for visualization.
  !****
  subroutine unnormalize (cell, dir, atoms, status)
    type(supercell), intent(inout) :: cell
    integer, intent(in) :: dir, atoms(:)
    integer, intent(out), optional :: status
    integer :: i
    real(wp) :: cc

    if (dir < 1 .or. dir > 3) then
       call error_msg (1, "Direction to unnormalize in must be 1, 2 or 3.", status)
       return
    end if
    if (cell%cartesian) then
       call cartesian2Direct(cell)
    end if
    if (.not. cell%relative) then
       call act2Rel(cell)
    end if
    do i = 1, size(atoms)
        cc = cell%atomCoords(dir, atoms(i))
        if (cc < 0.5_wp) then
            cell%atomCoords(dir, atoms(i)) = cc + 1.0_wp
        else
            cell%atomCoords(dir, atoms(i)) = cc - 1.0_wp
        end if
    end do
  end subroutine unnormalize


  !****f* supercell_modify/centercell
  ! PURPOSE
  ! Center a supercell around a point. Currently only support 0 and center.
  ! NOTES
  ! The idea is to for each atom whose coordinates are more than 0.5
  ! away from the chose point among some axis, test if moving the atom to
  ! its periodic image location will shorten the distance, and if so,
  ! update the coordinate.
  !****
  subroutine centercell (infile, point, outfile)
    character(len=*), intent(in) :: infile, outfile, point
    type(supercell) :: cell
    real(wp) :: pcoord(3), dist, tcoord(3)
    integer :: ii, jj

    call read_POSCAR (cell, infile)
    call cartesian2Direct (cell)
    if (point == "c") then
       pcoord = (/ 0.5_wp, 0.5_wp, 0.5_wp /)
    else
       pcoord = (/ 0.0_wp, 0.0_wp, 0.0_wp /)
    end if
    do ii = 1, size (cell%atoms)
       tcoord = cell%atomCoords(:,ii) - pcoord
       dist = sqrt (dot_product (tcoord, tcoord))
       do jj = 1, 3
          if (cell%atomCoords(jj,ii) > 0.5_wp) then
             tcoord = cell%atomCoords(:,ii)
             tcoord(jj) = tcoord(jj) - 1.0_wp
             if (sqrt (dot_product (tcoord, tcoord)) < dist) then
                cell%atomCoords(:,ii) = tcoord
             end if
          end if
       end do
    end do
    call write_POSCAR (cell, outfile)
  end subroutine centercell


  !****f* supercell_modify/interpolate_POSCAR
  ! PURPOSE
  !   Interpolate coordinates between two POSCAR files.
  !****
  subroutine interpolate_POSCAR(infile1, infile2, fraction, outfile)
    character(len=*), intent(in) :: infile1, infile2, outfile
    type(supercell) :: cell1, cell2
    real(wp), intent(in) :: fraction
    call read_POSCAR(cell1, infile1)
    call read_POSCAR(cell2, infile2)
    if (.not. checkCells(cell1, cell2)) then
       print *, 'Error: Cells are not consistent!'
       return
    end if
    call interpolate(cell1, cell2, fraction)
    call write_POSCAR(cell1, outfile)
  end subroutine interpolate_POSCAR


  !****f* supercell_core/interpolate
  ! PURPOSE
  ! Interpolate atom coordinates between two supercells.
  ! The result is returned in the first cell, whose coordinates
  ! are replaced with the interpolated ones.
  !****
  subroutine interpolate(cell1, cell2, fraction)
    type(supercell), intent(inout) :: cell1
    type(supercell), intent(inout) :: cell2
    real(wp), intent(in) :: fraction
    real(wp), dimension(3, size(cell1%atoms)) :: dist

    call cartesian2Direct(cell1)
    call cartesian2Direct(cell2)
    call act2rel(cell1)
    call act2rel(cell2)

    !      call normalize(cell1)
    !      call normalize(cell2)
    dist = cell2%atomCoords - cell1%atomCoords

    ! Take periodic bc into account
    where (dist > 0.5_wp)
       dist = dist - 1.0_wp
    end where
    where (dist < -0.5_wp)
       dist = dist + 1.0_wp
    end where

    cell1%atomCoords = cell1%atomCoords + fraction * dist
    call normalize(cell1)
    call check_nndist(cell1)
  end subroutine interpolate


  !****f* supercell_modify/removeAtoms_POSCAR
  ! PURPOSE
  !   Remove atoms from a POSCAR file.
  !****
  subroutine removeAtoms_POSCAR(infile, outfile, atomnums)
    character(len=*), intent(in) :: infile, outfile
    integer, dimension(:), intent(in) :: atomnums
    type(supercell) :: cell
    call read_POSCAR(cell, infile)
    print *, 'Removing ', size(atomnums), ' atoms.'
    call removeAtoms(cell, atomnums)
    call write_POSCAR(cell, outfile)
  end subroutine removeAtoms_POSCAR


  !****f* supercell_core/removeAtoms
  ! PURPOSE
  ! Remove atoms specified by their indexes starting from 1.
  !****
  subroutine removeAtoms(cell, atomnums)
    type(supercell), intent(inout) :: cell
    integer, dimension(:), intent(in) :: atomnums
    real(wp), dimension(3, size(cell%atoms) - size(atomnums)) :: atomCoords
    ! Portland pgf90 5.2-4 chokes on the line below, allocate manually
    ! instead.
    !type(atom), dimension(size(cell%atoms) - size(atomnums)) :: atoms
    type(atom), dimension(:), allocatable :: atoms
    logical, dimension(size(cell%atoms)) :: mask
    integer :: i, j

    !pgf90 fix
    allocate (atoms(size(cell%atoms) - size(atomnums)))
    mask = .TRUE.
    do i = 1, size(atomnums)
       mask(atomnums(i)) = .FALSE.
    end do
    j = 1
    do i = 1, size(cell%atoms)
       if (mask(i)) then
          atoms(j) = cell%atoms(i)
          atomCoords(:,j) = cell%atomCoords(:,i)
          j = j+1
       end if
    end do
    i = size(cell%atoms)-size(atomnums) ! atoms in new cell
    deallocate(cell%atoms, cell%atomCoords)
    allocate(cell%atoms(i), cell%atomCoords(3,i))
    cell%atoms = atoms
    cell%atomCoords = atomCoords
  end subroutine removeAtoms


  !****f* supercell_modify/lock_atoms
  ! PURPOSE
  ! Switch on selective dynamics and lock some atoms.
  !****
  subroutine lock_atoms (infile, outfile, atomnums)
    character(len=*), intent(in) :: infile, outfile
    integer, intent(in) :: atomnums(:)
    type (supercell) :: cell
    integer :: i
    call read_POSCAR (cell, infile)
    cell%selective = .true.
    do i = 1, size(atomnums)
       cell%atoms(atomnums(i))%selective = .false.
    end do
    call write_POSCAR (cell, outfile)
  end subroutine lock_atoms


end module supercell_modify
