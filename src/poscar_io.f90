!****h* vasputil/poscar_io
! NAME
!   poscar_io
! COPYRIGHT
!  Copyright (c) 2004, 2005, 2006, 2007, 2008 Janne Blomqvist

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
!   This module defines I/O helper procedures for 
!   reading and writing supercells in POSCAR format.
!****
module poscar_io
  use conf
  use supercell_core
  use supercell_utils

  implicit none

contains


  !****f* poscar_io/read_POSCAR
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


  !****f* poscar_io/write_POSCAR
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
       !         do i = 1, natoms
       !            print *, cell%atoms(i)%symbol
       !         end do
       print *, 'species: ', species
       do i = 1, size(symbols)
          print *, 'symbols: ', symbols(i)
       end do

    end if

    open(unit=pos_iu, file=outfile, form='formatted', &
         access='sequential', action='write', status='replace')
    write(pos_iu, '(100(A, 1X))') symbols
    write(pos_iu, *) cell%lattice%a
    write(pos_iu, '(3(4X, F18.15))') cell%lattice%t
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
          write(pos_iu, '(3(2X, F18.15), 3(2X, L2))') &
               cell%atomCoords(:,i), &
               cell%atoms(i)%selective
       end do
    else
       write(pos_iu, '(3(2X, F18.15))') (cell%atomCoords(:, i), i=1, size(cell%atoms))
    end if
    close(pos_iu)
    print *, trim(outfile), ' written successfully.'
  end subroutine write_POSCAR


end module poscar_io
