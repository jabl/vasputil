!****h* vasputil/sc_file_convert
! NAME
!   sc_file_convert
! COPYRIGHT
!  Copyright (c) 2004, 2005, 2006 Janne Blomqvist

!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.

!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

! PURPOSE

! This module defines procedures for converting between different
! file formats, and different variants of same formats (e.g. cartesian 
! vs. direct coordinates).

!****
module sc_file_convert
  use conf
  use poscar_io
  use supercell_utils

  implicit none

contains

  !****f* sc_file_convert/POSCAR2xyz
  ! NAME
  !   POSCAR2xyz
  ! PURPOSE
  !   Convert a POSCAR file to .xyz format.
  !****
  subroutine POSCAR2xyz(infile, outfile)
    character(len=*), intent(in) :: infile, outfile
    type(supercell) :: cell
    call read_POSCAR(cell, infile)
    !      call normalize(cell)
    call write_xyz(cell, outfile)
  end subroutine POSCAR2xyz


  !****f* supercell_io/direct2cartesian_poscar
  ! PURPOSE
  ! Convert POSCAR from direct to cartesian coordinates.
  !****
  subroutine direct2cartesian_poscar (infile, outfile)
    character(len=*), intent(in) :: infile, outfile
    type(supercell) :: cell
    call read_POSCAR (cell, infile)
    call direct2Cartesian (cell)
    call write_POSCAR (cell, outfile)
  end subroutine direct2cartesian_poscar


  !****f* sc_file_convert/kspace2xyz
  ! PURPOSE
  !   Convert a POSCAR file to .xyz format, in reciprocal coordinates.
  !****
  subroutine kspace2xyz(infile, outfile)
    character(len=*), intent(in) :: infile, outfile
    type(supercell) :: cell
    call read_POSCAR(cell, infile)
    call direct2Cartesian(cell)
    !      call rel2Act(cell)
    call direct2recip(cell%lattice, cell%atomCoords)
    call write_xyz(cell, outfile)
  end subroutine kspace2xyz


  !****f* sc_file_convert/xyz2poscar
  ! PURPOSE
  ! Convert a xyz file to POSCAR format.
  !****
  subroutine xyz2poscar (infile, outfile)
    character(len=*), intent(in) :: infile, outfile
    type(supercell) :: cell
    call read_xyz (cell, infile)
    call write_POSCAR (cell, outfile)
  end subroutine xyz2poscar


  !****f* sc_file_convert/dumpcoords
  ! PURPOSE
  ! Dump cartesian coordinates to stdout.
  !****
  subroutine dumpcoords (infile)
    character(len=*), intent(in) :: infile
    type(supercell) :: cell
    call read_POSCAR (cell, infile)
    call direct2Cartesian (cell)
    call rel2Act (cell)
    print '(3F18.14)', cell%atomCoords
  end subroutine dumpcoords


  !****f* sc_file_convert/importcoords
  ! PURPOSE
  ! Import cartesian coordinates from stdin.
  !****
  subroutine importcoords (infile)
    character(len=*), intent(in) :: infile
    type(supercell) :: cell
    integer :: i, k
    call read_POSCAR (cell, infile)
    call direct2Cartesian (cell)
    call rel2Act (cell)
    !      read (*, '(3F13.8)') (cell%atomCoords (:, i), i=1,3)
    !      read (*, *) (cell%atomCoords (:, i), i=1,3)
    do k = 1, size (cell%atoms)
       read (*, *) &
            (cell%atomCoords(i, k), i=1,3)
    end do
    call write_POSCAR (cell, infile)
  end subroutine importcoords


  !****f* sc_file_convert/dumpatomsase
  ! PURPOSE
  ! Dump cartesian coordinates and atom types to stdout
  ! in a format that can be understood by Campos ASE
  ! (i.e. as python code).
  !****
  subroutine dumpatomsase (infile)
    character(len=*), intent(in) :: infile
    type(supercell) :: cell
    integer :: i, j
    logical :: writecomma
    call read_POSCAR (cell, infile)
    call direct2Cartesian (cell)
    call rel2Act (cell)
    write (*,'(A)') 'from ASE import Atom, ListOfAtoms'
    write (*,'(A)') 'atoms = ListOfAtoms(['
    writecomma = .false.
    do i = 1, size (cell%atoms)
       if (writecomma) then
          write (*,*) ','
       end if
       write (*,'(A,A,A,3(F13.8,A),A)', advance="no") 'Atom("', &
            trim (cell%atoms(i)%symbol), '", (', &
            cell%atomcoords(1,i), ', ', cell%atomcoords(2,i), ', ', &
            cell%atomcoords(3,i), '))'
       writecomma = .true.
    end do
    write (*,*) '], cell=('
    writecomma = .false.
    do i = 1, 3
       if (writecomma) then
          write (*,*) ','
       end if
       write (*,'(A)', advance="no") '('
       writecomma = .false.
       do j = 1, 3
          if (writecomma) then
             write (*,'(A)', advance="no") ', '
          end if
          write (*,'(F13.8)', advance="no") cell%lattice%t(j,i)*cell%lattice%a
          writecomma = .true.
       end do
       write (*,*) ') '
       writecomma = .true.
    end do
    write (*,*) '), periodic=1)'
  end subroutine dumpatomsase


  !****f* sc_file_convert/read_xyz
  ! PURPOSE
  ! Read a .xyz file.
  !****
  subroutine read_xyz (cell, infile)
    character(len=*), intent(in) :: infile
    type(supercell), intent(out) :: cell
    integer, parameter :: xyz_iu = 12
    integer :: natoms, i
    character(len=5) :: symbol

    open (unit=xyz_iu, file=infile, form='formatted', &
         access='sequential', action='read', status='old')
    read (xyz_iu, *) natoms
    read (xyz_iu, '(A)') ! Comment line
    call init_cell (cell, natoms)
    cell%cartesian = .true.
    cell%relative = .false.
    do i = 1, natoms
       read (xyz_iu, '((A5), 3(2X, F12.8))') symbol, cell%atomCoords(:,i)
       cell%atoms(i)%symbol = adjustl (symbol)
    end do
    close (xyz_iu)
    call read_lattice (cell%lattice)
    call sort_cell (cell)
  end subroutine read_xyz


  !****f* sc_file_convert/write_xyz
  ! PURPOSE
  !   Write a .xyz file using the supplied supercell.
  !****
  subroutine write_xyz(cell, outfile)
    character(len=*), intent(in) :: outfile
    type(supercell), intent(inout) :: cell
    integer, parameter :: xyz_iu = 12 ! IO unit for .xyz file
    integer :: i

    call direct2Cartesian(cell)
    call rel2Act(cell)
    open(unit=xyz_iu, file=outfile, form='formatted', &
         access='sequential', action='write', status='replace')
    write(xyz_iu, '(I0)') size(cell%atoms)
    write(xyz_iu, '(A)') 'Generated by vasputil, (C) 2004, 2005, 2006 Janne Blomqvist'
    do i = 1, size(cell%atoms)
       write(xyz_iu, '((A), 3(2X, F12.8))') cell%atoms(i)%symbol, cell%atomCoords(:,i)
    end do
    close(xyz_iu)
    print *, trim(outfile), ' written successfully.'
  end subroutine write_xyz


end module sc_file_convert
