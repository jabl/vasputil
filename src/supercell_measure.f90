!****h* vasputil/supercell_measure
! NAME
!   supercell_measure
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
module supercell_measure
  use conf
  use geometry
  use lattice
  use poscar_io
  use supercell_utils

  implicit none

contains


  !****f* supercell_measure/atomdist
  ! PURPOSE
  ! Calculate the distance between two atoms in a supercell.
  ! This function takes into account periodic boundary conditions.
  !****
  function atomdist(cell1, atom1, cell2, atom2, errmsg)
    type(supercell), intent(inout), target :: cell1
    type(supercell), intent(inout), optional, target :: cell2
    integer :: atom1, atom2
    character(len=*), intent(out), optional :: errmsg
    real(wp) :: atomdist
    type(supercell), pointer :: cello
    real(wp), dimension(3) :: distVect

    if (present (errmsg)) then
       errmsg = ""
    end if
    if (present (cell2)) then
       cello => cell2
       if (.not. checkCells(cell1, cell2, errmsg)) then
          return
       end if
       call cartesian2Direct (cello)
       call act2Rel (cello)
    else
       cello => cell1
    end if

    call cartesian2Direct (cell1)
    call act2Rel (cell1)
    ! Sanity check
    if (atom1 > size (cell1%atoms) .or. atom2 > size (cello%atoms)) then
       if (present (errmsg)) then
          errmsg = "Atom index larger than the number of atoms in the cell"
       end if
       return
    end if

    distVect = cell1%atomCoords(:, atom1) - cello%atomCoords(:, atom2)
    print *, 'cell1 atom: ', cell1%atomCoords(:, atom1)
    print *, 'cell2 atom: ', cello%atomCoords(:, atom2)
    print *, 'distvect is: ', distvect
    distVect = distVect - anint (distVect)
    distVect = matmul (cell1%lattice%t, distVect)
    atomdist = sqrt (dot_product (distVect, distVect)) * cell1%lattice%a
  end function atomdist


  !****f* supercell_measure/planetopoint
  ! PURPOSE
  ! Calculates the average distance from a plane to a set of atoms
  !****
  function planetopoint(cell, planeatom1, planeatom2, planeatom3, atoms)
    type(supercell), intent(inout) :: cell
    integer, intent(in) :: planeatom1, planeatom2, planeatom3
    integer, dimension(:), intent(in) :: atoms
    type(plane3) :: plane
    real(wp) :: planetopoint, tot_dist
    integer :: i

    !      call normalize(cell)
    call direct2Cartesian(cell)
    call rel2Act(cell)
    call new(plane, cell%atomCoords(:, planeatom1), &
         cell%atomCoords(:, planeatom2), cell%atomCoords(:, planeatom3))
    tot_dist = 0._wp
    do i = 1, size(atoms)
       tot_dist = tot_dist + distanceFrom(plane, cell%atomCoords(:, atoms(i)))
    end do
    planetopoint = tot_dist/size(atoms)
  end function planetopoint


  !****f* supercell_measure/planetoatom
  ! PURPOSE
  !   Calculate the average distance from a plane to a set of atoms.
  !****
  subroutine planetoatom(infile, pa1, pa2, pa3, atoms)
    character(len=*), intent(in) :: infile
    type(supercell) :: cell
    integer, intent(in) :: pa1, pa2, pa3
    integer, dimension(:), intent(in) :: atoms
    call read_POSCAR(cell, infile)
    print *, '(Average) distance from plane to atom(s) is: ', &
         planetopoint(cell, pa1, pa2, pa3, atoms)
  end subroutine planetoatom


  !****f* supercell_measure/planetolayer
  ! PURPOSE
  !   Calculate the distance from a plane to a layer.
  !****
  subroutine planetolayer(infile, pa1, pa2, pa3, lheight, tol, bulk)
    character(len=*), intent(in) :: infile
    type(supercell) :: cell
    integer, intent(in) :: pa1, pa2, pa3
    real(wp), intent(in) :: lheight
    real(wp), intent(in), optional :: tol, bulk
    real(wp) :: t
    logical, dimension(:), allocatable :: layermask
    integer, dimension(:), allocatable :: atoms
    integer :: i, j
    call read_POSCAR(cell, infile)
    if (present(tol)) then
       t = tol
    else
       t = 0.5_wp
    end if
    allocate(layermask(size(cell%atoms)))
    layermask = .FALSE.
    call direct2Cartesian(cell)
    call rel2Act(cell)
    !      print *, cell%atomCoords(:, 3)
    where (abs(cell%atomCoords(3, :) - lheight) < t)
       layermask = .TRUE.
    end where
    !      print *, layermask
    allocate(atoms(count(layermask)))
    j = 1
    do i = 1, size(layermask)
       if (layermask(i)) then
          atoms(j) = i
          j = j+1
       end if
    end do
    t = planetopoint( cell, pa1, pa2, pa3, atoms)
    print *, 'Total of ', size(atoms), ' atoms included in layer.'
    print *, 'Atom numbers included in layer are: ', atoms
    print *, 'Average distance from plane to layer is ', t
    if (present(bulk)) then
       print '(A, F7.2, A)', 'Relaxation compared to bulk is: ', (t/bulk - 1.0_wp)*100.0_wp, ' %.'
    end if
  end subroutine planetolayer


  !****f* supercell_measure/atomsDistance
  ! PURPOSE
  !   Calculate the distance between two atoms.
  !****
  subroutine atomsDistance(infile1, atom1, infile2, atom2)
    character(len=*), intent(in) :: infile1
    character(len=*), intent(in), optional :: infile2
    integer, intent(in) :: atom1, atom2
    type(supercell) :: cell1, cell2
    real(wp) :: dist
    character(len=132) :: errmesg = ''
    call read_POSCAR(cell1, infile1)
    if (present (infile2)) then
       if (infile1 /= infile2) then
          call read_POSCAR(cell2, infile2)
          dist = atomdist (cell1, atom1, cell2=cell2, atom2=atom2, errmsg=errmesg) 
       else
          dist = atomdist (cell1, atom1, atom2=atom2, errmsg=errmesg)
       end if
    else
       dist = atomdist (cell1, atom1=atom1, atom2=atom2, errmsg=errmesg)
    end if
    if (errmesg /= "") then
       call error_stop (errmesg)
    else
       print *, 'Distance between atoms is: ', dist
    end if
  end subroutine atomsDistance


  !****f* supercell_measure/atomsMoved
  ! FUNCTION
  ! Compare two POSCAR files and print out a list of atom indexes for those
  ! atoms that have moved more than the specified distance.
  !****
  subroutine atomsMoved(infile1, infile2, tol)
    character(len=*), intent(in) :: infile1, infile2
    real(wp), intent(in), optional :: tol
    real(wp) :: t, dist
    type(supercell) :: cell1, cell2
    integer :: i

    call read_POSCAR(cell1, infile1)
    call read_POSCAR(cell2, infile2)
    if (.not. checkCells(cell1, cell2)) then
       print *, 'Error: Cells not consistent!'
       return
    end if
    if (present(tol)) then
       t = tol
    else
       t = 0.1_wp
    end if
    call direct2Cartesian(cell1)
    call direct2Cartesian(cell2)
    call rel2Act(cell1)
    call rel2Act(cell2)
    print *, 'Atoms that have moved more than ', t, ' Å:'
    do i = 1, size(cell1%atoms)
       dist = norm(cell1%atomCoords(:,i) - cell2%atomCoords(:,i))
       if (dist > t) then
          print *, i, ': ', dist, ' Å'
       end if
    end do
  end subroutine atomsMoved



  !****f* supercell_measure/check_nn_POSCAR
  ! FUNCTION
  ! Check nearest neighbor distances between atoms in a supercell.
  !****
  subroutine check_nn_POSCAR(infile, tol)
    character(len=*), intent(in) :: infile
    real(wp), intent(in), optional :: tol
    type(supercell) :: cell

    call read_POSCAR(cell, infile)
    if (present(tol)) then
       call check_nndist(cell, tol)
    else
       call check_nndist(cell)
    end if
  end subroutine check_nn_POSCAR


end module supercell_measure
