!****h* vasputil/supercell_io
! NAME
!   supercell_io
! COPYRIGHT
!  Copyright (c) 2004, 2005 Janne Blomqvist

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
!   This module defines I/O helper procedures for 
!   reading and writing supercells.
!****
module supercell_io
  use conf
  use supercell_core

  implicit none

  contains

    !****f* supercell_io/POSCAR2xyz
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


    !****f* supercell_io/kspace2xyz
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


    !****f* supercell_io/xyz2poscar
    ! PURPOSE
    ! Convert a xyz file to POSCAR format.
    !****
    subroutine xyz2poscar (infile, outfile)
      character(len=*), intent(in) :: infile, outfile
      type(supercell) :: cell
      call read_xyz (cell, infile)
      call write_POSCAR (cell, outfile)
    end subroutine xyz2poscar


    !****f* supercell_io/planetoatom
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


    !****f* supercell_io/planetolayer
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

    
    !****f* supercell_io/atomsDistance
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


    !****f* supercell_io/normalize_POSCAR
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


    !****f* supercell_io/unnormalize_POSCAR
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


    !****f* supercell_io/dumpcoords
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


    !****f* supercell_io/importcoords
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


    !****f* supercell_io/interpolate_POSCAR
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


    !****f* supercell_io/atomsMoved
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

    
    !****f* supercell_io/removeAtoms_POSCAR
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


    !****f* supercell_io/check_nn_POSCAR
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


    !****f* supercell_io/lock_atoms
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


    !****f* supercell_io/sc_generator_io
    ! NAME
    !   sc_generator_io -- Generate a supercell
    !****
    subroutine sc_generator_io (infile, infile2, outfile)
      character(len=*), intent(in) :: infile, outfile
      character(len=*), optional, intent(in) :: infile2
      type(supercell) :: cell, scell
      integer :: stat

      if (present (infile2)) then
         call read_POSCAR (scell, infile2, stat)
         !call write_POSCAR (scell, "test")
      else
         call read_lattice (scell%lattice)
      end if
      ! Read the basis of the cell that is to be replicated in the 
      ! new supercell.
      call read_POSCAR (cell, infile)
      ! Fill out the new supercell
      call generate_supercell (cell, scell, error_stop)
      ! Write the new supercell
      call cartesian2Direct (scell)
      call act2rel (scell)
      call write_POSCAR (scell, outfile)
    end subroutine sc_generator_io

      
    !****f* supercell_io/read_POSCAR
    ! PURPOSE
    ! Read a POSCAR file into a derived type.
    !****
    subroutine read_POSCAR(cell, infile, status)
      character(len=*), intent(in) :: infile
      type(supercell), intent(out) :: cell
      integer, intent(out), optional :: status
      integer, parameter :: pos_iu = 10, & ! IO unit for POSCAR
           pot_iu = 11 ! IO unit for POTCAR
      character(len=132) :: first, line
      integer :: i, j, spec_counter, k, lower, upper, stat
      logical :: spc, & ! Last before current character in string is a space?
           pot_ex 
      integer, dimension(:), allocatable :: species

      cell%relative = .TRUE.
      open(unit=pos_iu, file=infile, form='formatted', &
           access='sequential', action='read', status='old')
      read(pos_iu, '(A)') first ! First line, a comment?
      read (pos_iu, '(F16.10)') cell%lattice%a
!      read(pos_iu, '(3(5X, F16.1))') (cell%latticeVectors(j,:), j=1,3)
      read(pos_iu, *) (cell%lattice%t(:, j), j=1,3)
      if (is_lin_dependent (cell%lattice%t)) then
         call error_stop ("ERROR: Lattice vectors are not linearly independent!")
      end if
      !print *, 'Lattice(2,1): ', cell%lattice%t(2,1)
      spec_counter = 0 ! how many atom species are there?
      read(pos_iu, '(A)') line
      line = adjustl(line)
!      print *, len_trim(line)
!      line = trim(line)
      spc = .TRUE.
      do i = 1,len_trim(line)
         if ( line(i:i) /= ' ' .and. spc ) then
            spec_counter = spec_counter + 1
            spc = .FALSE.
         else if ( line(i:i) == ' ') then
            spc = .TRUE.
         end if
      end do
      allocate(species(spec_counter))
!      species = 0
!      print *, spec_counter, ' atom species found.'
      read(line, *) species
!      print *, 'Atoms of each species: ', species
      read(pos_iu, '(A)') line
      if (line(1:1) == 'S' .or. line(1:1) == 's') then ! Selective dynamics
         cell%selective = .TRUE.
         read(pos_iu, '(A)') line ! Cartesian
      else
         cell%selective = .FALSE.
      end if
      if (line(1:1) == 'C' .or. line(1:1) == 'c' .or. &
           line(1:1) == 'K' .or. line(1:1) == 'k') then
         cell%cartesian = .TRUE.
      else
         cell%cartesian = .FALSE.
      end if
      ! Allocate memory for the atom coordinates and info
      j = sum (species) ! How many atoms are there in total
      if (debug) then
         print *, 'Total of ', j, ' atoms of ', &
              spec_counter, ' species found: ', species
      end if
      !allocate( cell%atoms(j), cell%atomCoords(3, j))
      call init_cell (cell, j)
      ! Read in the atom coordinates
      do k = 1, j
         if (cell%selective) then
!            read(pos_iu, '(3(2X, F16.10, 2X), 3(2X, L2))') &
            read (pos_iu, *) &
                 (cell%atomCoords(i, k), i=1,3), &
                 (cell%atoms(k)%selective(i), i=1,3)
         else
!            read(pos_iu, '(3(2X, F16.10, 2X))') &
            read (pos_iu, *) &
                 (cell%atomCoords(i, k), i=1,3)
         end if
      end do
      close(pos_iu)
!!$      do i=1,sum(species)
!!$         print *, cell%atomCoords(i,:)
!!$         print *, cell%atoms(i)
!!$      end do

      ! Try to determine the atomic symbols of the atoms in the cell.
      ! We have three ways to do this, we try in this order:
      ! 1. Check if POTCAR exists and take symbols from it.
      ! 2. Check if the first line of POSCAR contains them.
      ! 3. Ask the user to supply them. 

      ! See if POTCAR exists
      inquire (file='POTCAR', exist= pot_ex)
      if (pot_ex) then
         open (unit=pot_iu, file='POTCAR', form='formatted', &
              access='sequential', action='read', status='old')
         j = 1
         k = 0 !offset
         do while (j <= size (species))
            read(pot_iu, '(A)') line
            if ( index (line, 'VRHFIN') == 0) then
               cycle
            end if
            upper = index (line, ':') - 1
            lower = index (line, '=') + 1
            do i = 1, species(j)
               cell%atoms(i+k)%symbol = trim (adjustl (line(lower:upper)))
            end do
            k = k + species(j)
            j = j + 1
         end do
         close(pot_iu)

      else ! Check the first line
         call set_species (cell%atoms, species, first, stat)
         if (stat /= 0 .and. (.not. present (status))) then
            print *, "Couldn't read species from ", trim (infile), &
                 ", status code: ", stat
            ! Some error has occured, ask the user instead
            call read_species (cell%atoms, species)
         else
            if (present (status)) then
               status = stat
            end if
         end if
      end if

      call latt_init (cell%lattice)

      if (debug) then
         print *, trim (infile), ' read successfully.'
      end if

    contains

      
      !****f* read_POSCAR/eat_whitespace
      ! PURPOSE
      ! Increase the index until a non-whitespace character is encountered.
      !****
      subroutine eat_whitespace (string, cursor, status)
        character(len=*), intent(in) :: string
        integer, intent(inout) :: cursor
        integer, intent(out), optional :: status

        if (present (status)) then
           status = 0
        end if
        if (cursor >= len (string) .or. cursor < 0) then
           call error_msg (1, 'Cursor is not within string bounds', status)
        end if
        do while (string(cursor:cursor) == ' ')
           cursor = cursor + 1
           if (cursor > len (string)) then
              call error_msg (2, 'Encountered end of string', status)
              exit
           end if
        end do
      end subroutine eat_whitespace


      !****f* read_POSCAR/read_species
      ! PURPOSE
      ! Read the atom species from the terminal.
      !****
      subroutine read_species (atoms, species)
        type(atom), dimension(:), intent(inout) :: atoms
        integer, intent(in) :: species(:)
        character(len=132) :: symbols
        integer :: stat

        print *, 'Enter the atomic symbols of the ', size(species), &
             ' species in the supercell, &
             &separated by space:'
        read (*, '(A)') symbols
        call set_species (atoms, species, symbols, stat)
        if (stat /= 0) then
           call error_stop ('Error occured when trying to read species from terminal')
        end if
      end subroutine read_species


      !****f* read_POSCAR/set_species
      ! PURPOSE
      ! Set the atomic species from a string.
      !****
      subroutine set_species (atoms, species, string, status)
        type(atom), intent(inout) :: atoms(:)
        integer, intent(in) :: species(:)
        character(len=*), intent(in) :: string
        integer, intent(out) :: status
        integer :: lower, upper, i, j, k, stat, strlen

        status = 0
        lower = 1
        k = 0 ! offset index
        do i = 1, size(species)
           call eat_whitespace (string, lower, stat)
           !print *, 'lower after eat_whitespace: ', lower
           upper = index(string(lower:), ' ') + lower - 1
           strlen = upper - lower
           !print *, '$', string(lower:upper - 1), '$', strlen
           if (strlen < 1 .or. strlen > 2 .or. stat /= 0) then
              !print *, 'strlen: ', strlen, ' stat: ', stat
              !print *, string(lower:upper)
              status = 1
              return
           end if
           do j = 1, species(i)
              atoms(j+k)%symbol = trim(string(lower:upper - 1))
           end do
           k = k + species(i)
           lower = upper
        end do
      end subroutine set_species

    end subroutine read_POSCAR


    !****f* supercell_io/write_POSCAR
    ! PURPOSE
    !   Write a POSCAR file using the supplies supercell.
    !****
    subroutine write_POSCAR(cell, outfile)
      character(len=*), intent(in) :: outfile
      type(supercell), intent(inout) :: cell
      integer, parameter :: pos_iu = 13 ! IO unit for output file
      integer :: i, j, spec, natoms
      integer, dimension(:), allocatable :: species
      character(len=2), dimension(:), allocatable :: symbols

      call act2Rel(cell) ! Make sure coordinates are relative

!      do i = 1, size(cell%atoms)
!         print *, cell%atoms(i)%symbol
!      end do

      ! First determine how many atomic species we have
      ! This algorithm assumes that the species are ordered.
      natoms = size(cell%atoms)
      if (natoms <= 1) then
         spec = natoms
      else if (natoms > 1) then
         spec = 1
         do i = 2, natoms
            if (cell%atoms(i)%symbol /= cell%atoms(i-1)%symbol) then
               spec = spec + 1
            end if
         end do
      end if

!      print *, 'Allocating for ', spec, ' species.'
      allocate( species(spec), symbols(spec))
      species = 0
      if (natoms == 1) then
         species = 1
         symbols(1) = cell%atoms(1)%symbol
      else if (natoms > 1) then
         species = 1
         j = 1
         symbols(1) = cell%atoms(1)%symbol
         do i = 2, natoms
            if (cell%atoms(i)%symbol == cell%atoms(i-1)%symbol) then
               species(j) = species(j) + 1
            else
               j = j+1
               symbols(j) = cell%atoms(i)%symbol
            end if
         end do
      end if

      if (debug) then
         print *, 'Atoms: ', natoms
         do i = 1, natoms
            print *, cell%atoms(i)%symbol
         end do
         print *, 'species: ', species
         do i = 1, size(symbols)
            print *, 'symbols: ', symbols(i)
         end do

      end if

      open(unit=pos_iu, file=outfile, form='formatted', &
           access='sequential', action='write', status='replace')
      write(pos_iu, '(100(A, 1X))') symbols
      write(pos_iu, *) cell%lattice%a
      write(pos_iu, '(3(1X, F18.14))') cell%lattice%t
      write(pos_iu, *) species
      if (cell%selective) then
         write(pos_iu, '(A)') 'Selective dynamics'
      end if
      if (cell%cartesian) then
         write(pos_iu, '(A)') 'Cartesian'
      else
         write(pos_iu, '(A)') 'Direct'
      end if
      if (cell%selective) then
         do i = 1, size(cell%atoms)
            write(pos_iu, '(3(2X, F18.14, 2X), 3(2X, L2))') &
                 cell%atomCoords(:,i), &
                 cell%atoms(i)%selective
         end do
      else
         write(pos_iu, '(3(1X, F18.14))') (cell%atomCoords(:, i), i=1, size(cell%atoms))
      end if
      close(pos_iu)
      print *, trim(outfile), ' written successfully.'
    end subroutine write_POSCAR


    !****f* supercell_io/read_xyz
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

    !****f* supercell_io/write_xyz
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
      write(xyz_iu, '(A)') 'Generated by vasputil, (C) 2004, 2005 Janne Blomqvist'
      do i = 1, size(cell%atoms)
         write(xyz_iu, '((A), 3(2X, F12.8))') cell%atoms(i)%symbol, cell%atomCoords(:,i)
      end do
      close(xyz_iu)
      print *, trim(outfile), ' written successfully.'
    end subroutine write_xyz


end module supercell_io
