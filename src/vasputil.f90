!****h* /vasputil
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
! VASP utils main program. It works by checking the 1:st command line
! argument for the name of the command to run, and from that
! determining what to do.
!****
program vasputil
  use conf
  use f2kcli
  use supercell_measure
  use supercell_modify
  use sc_file_convert
  use supercell_generator

  implicit none

  character(len=80) :: command
  character(len=132), dimension(:), allocatable :: arg
  integer :: status, i, j, iarg
  integer, dimension(3) :: pta_arg
  real(wp) :: fraction, lheight, tol, bulk
  integer, dimension(:), allocatable :: atoms

  iarg = command_argument_count ()
  
  if (iarg == 0) then
     call print_usage ()
     stop
  end if

  ! Get the command name     
  call get_command_argument (number=1,  value=command, status=status)
  if ( status /= 0) then
     print *, 'Could not retrieve name of command, status: ', status
     stop
  end if

  ! Allocate space for the arguments. We need two extra spaces for non-present
  ! optional args.
  allocate(arg(iarg+1))

  do j = 2, iarg
     call get_command_argument (j, arg(j-1), status=status)
     if (status /= 0) then
        print *, 'Could not retrieve value of command argument ', j, &
             ', status: ', status
        stop
     end if
  end do

  i = iarg - 1

  select case (command)
  case ('poscar2xyz')
     select case (i)
     case default
        print *, 'Usage: poscar2xyz POSCAR [POSCAR.xyz]'
        print *, ''
        print *, 'Convert a POSCAR file to xyz format.'
        print *, 'POSCAR is the POSCAR file you want to convert.'
        print *, 'POSCAR.xyz is the filename of the .xyz file to be &
             &created. If omitted, it defaults to the name of the &
             &POSCAR file, plus the xyz extension.'
     case (1, 2)
        if ( i == 1) then
           write(arg(2), '(A)') trim (arg(1)) // '.xyz'
        end if
        call POSCAR2xyz (arg(1), arg(2))
     end select

  case ('direct2cartesian')
     select case (i)
     case default
        print *, 'Usage: direct2cartesian POSCAR [POSCAR.cart]'
        print *, ''
        print *, 'Convert a POSCAR file in direct coordinates to cartesian.'
        print *, 'POSCAR: The POSCAR file with the direct coordinates.'
        print *, 'POSCAR.cart: The name of the output file. If omitted, &
             &defaults to the name of the POSCAR file plus the .cart extension.'
     case (1, 2)
        if (i == 1) then
           write (arg(2), '(A)') trim (arg(1)) // '.cart'
        end if
        call direct2cartesian_poscar (arg(1), arg(2))
     end select

  case ('xyz2poscar')
     select case (i)
     case default
        print *, 'Usage: xyz2poscar XYZfile POSCARfile'
        print *, ''
        print *, 'Convert a XYZ file to POSCAR format.'
        print *, 'XYZfile: The name of the .xyz format file.'
        print *, 'POSCARfile: The name of the output file.'
     case (2)
        call xyz2poscar (arg(1), arg(2))
     end select

  case ('atomsmoved')
     select case (i)
     case default
        print *, 'Usage: atomsmoved POSCAR1 POSCAR2 [tol]'
        print *, ''
        print *, 'Compare two POSCAR files to see which atoms have moved.'
        print *, 'POSCAR[1,2]: The file names of the files to compare.'
        print *, 'tol: Optional argument specifying the tolerance. ', &
             'If omitted, defaults to 0.1 Å'
     case (2, 3)
        if (i==2) then
           call atomsMoved (arg(1), arg(2))
        else
           read(arg(3), *) tol
           call atomsMoved (arg(1), arg(2), tol)
        end if
     end select

  case ('plane2atom')
     select case (i)
     case (:4)
        print *, 'Usage: plane2atom POSCAR planeatom1 planetatom2 planeatom3 atom1-N'
        print *, ''
        print *, 'Calculates the (average) distance between a plane and a set of atoms.'
        print *, 'POSCAR: the file containing the coordinates.'
        print *, 'planeatomX: The atoms defining the plane.'
        print *, 'atom1-N: The atom(s) to calculate the (average) distance to.'
     case default
        do j = 2,4
           read(arg(j), '(I3)') pta_arg(j-1)
        end do
        allocate(atoms(i-4))
        do j = 5, i
           read(arg(j), *) atoms(j-4)
        end do
        call planetoatom (arg(1), pta_arg(1), pta_arg(2), pta_arg(3), atoms)
     end select

  case ('plane2layer')
     select case (i)
     case (:4, 8:)
        print *, 'Usage: plane2layer POSCAR ptom1 patom2 patom3 &
             &layerheight [tol] [bulk]'
        print *, ''
        print *, 'Calculates the distance between a plane and a set of atoms &
             &determined by the value of their z coordinate.'
        print *, ''
        print *, 'POSCAR: the file containing the coordinates.'
        print *, ''
        print *, 'patomX: The atoms defining the plane.'
        print *, ''
        print *, 'layerheight: The approximate height (z coordinate) of the &
             &layer to calculate the distance to.'
        print *, ''
        print *, 'tol: Optional argument specifying the tolerance for the &
             &layerheight argument. If omitted, defaults to 0.5.'
        print *, ''
        print *, 'bulk: Optional argument specifying the bulk interlayer ', &
             'distance. If specified, the program prints out the layer ', &
             'relaxation in percent. The tol argument must also be present ', &
             'for this to be calculated.'
     case default
        do j = 2, 4
           read(arg(j), *) pta_arg(j-1)
        end do
        read(arg(5), *) lheight
        if (i == 6) then
           read(arg(6), *) tol
           call planetolayer (arg(1), pta_arg(1), pta_arg(2), pta_arg(3), lheight, tol)
        else if (i == 7) then
           read(arg(6), *) tol
           read(arg(7), *) bulk
           call planetolayer (arg(1), pta_arg(1), pta_arg(2), pta_arg(3), lheight, tol, bulk)
        else
           call planetolayer (arg(1), pta_arg(1), pta_arg(2), pta_arg(3), lheight)
        end if
     end select

  case ('atomsdistance')
     select case (i)
     case default
        print *, 'Usage: atomsdistance POSCAR1 atom1 [POSCAR2] atom2'
        print *, ''
        print *, 'Calculates the distance between two atoms.'
        print *, 'POSCAR1: The file containing the coordinates.'
        print *, 'POSCAR2: The optional file containing the coordinates ', &
             'of the second atom. If omitted, the coordinates are taken ', &
             'from the first file.'
        print *, 'atomX: the numbers of the atoms.'
     case (3, 4)
        read(arg(2), '(I3)') pta_arg(1)
        if (i == 3) then
           read(arg(3), '(I3)') pta_arg(2)
           arg(3) = arg(1)
        else
           read(arg(4), '(I3)') pta_arg(2)
        end if
        call atomsDistance (arg(1), pta_arg(1), arg(3), pta_arg(2))
     end select

  case ('interpolate')
     select case (i)
     case default
        print *, 'Usage: interpolate POSCAR1 POSCAR2 fraction POSCAR.out'
        print *, ''
        print *, 'Interpolates the coordinates between two POSCAR files.'
        print *, 'POSCARX: The two files to interpolate between.'
        print *, 'fraction: The fraction [0, 1] to interpolate with.'
        print *, 'POSCAR.out: The output file.'
     case (4)
        read(arg(3), *) fraction
        call interpolate_POSCAR (arg(1), arg(2), fraction, arg(4))
     end select


  case ('unnormalize')
     select case (i)
     case (:3)
        print *, 'Usage: unnormalize POSCAR POSCAR.out dir atom1-N'
        print *, ''
        print *, 'Unnormalize a POSCAR file, e.g. for visualization.'
        print *, 'POSCAR: The input file.'
        print *, 'POSCAR.out: The output file.'
        print *, 'dir: The lattice vector in which direction the ', &
             'coordinates should be changed.'
     case default
        read(arg(3), *) status  ! reuse status variable for dir
        allocate(atoms(i-3))
        do j = 4, i
           read(arg(j), *) atoms(j-3)
        end do
        call unnormalize_POSCAR (arg(1), status, atoms, arg(2))
     end select

  case ('normalize')
     select case (i)
     case default
        print *, 'Usage: normalize POSCAR [POSCAR.out]'
        print *, ''
        print *, 'Normalizes a POSCAR file.'
        print *, 'POSCAR: The input file.'
        print *, 'POSCAR.out: The output file. If omitted, defaults to ', &
             'the name of the input file plus the .out extension.'
     case (1, 2)
        if (i == 1) then
           write(arg(2), '(A)') trim(arg(1)) // '.out'
        end if
        call normalize_POSCAR (arg(1), arg(2))
     end select

  case ('removeatoms')
     select case (i)
     case (:2)
        print *, 'Usage: removeatoms POSCAR.in POSCAR.out atom[1-N]'
        print *, ''
        print *, 'Remove atoms from a POSCAR file.'
        print *, 'POSCAR.in: The input file.'
        print *, 'POSCAR.out: The output file.'
        print *, 'atom[1-N]: The indexes of the atoms to remove.'
     case default
        allocate(atoms(i-2))
        do j = 3, i
           read(arg(j), *) atoms(j-2)
        end do
        call removeAtoms_POSCAR (arg(1), arg(2), atoms)
     end select

  case ('lockatoms')
     select case (i)
     case (:2)
        print *, 'Usage: lockatoms POSCAR.in POSCAR.out atom[1-N]'
        print *, ''
        print *, 'Lock atom coordinates in a POSCAR file.'
        print *, 'POSCAR.in: The input file.'
        print *, 'POSCAR.out: The output file.'
        print *, 'atom[1-N]: The indexes of the atoms to lock.'
     case default
        allocate (atoms(i-2))
        do j = 3, i
           read (arg(j), *) atoms(j-2)
        end do
        call lock_atoms (arg(1), arg(2), atoms)
     end select

  case ('kspace2xyz')
     select case (i)
     case default
        print *, 'Usage: kspace2xyz POSCAR [POSCAR.k.xyz]'
        print *, ''
        print *, 'Convert a POSCAR file to xyz format with reciprocal &
             &coordinates.'
        print *, 'POSCAR is the POSCAR file you want to convert.'
        print *, 'POSCAR.k.xyz is the filename of the .xyz file to be'
        print *, 'created. If ommitted, it defaults to the name of the'
        print *, 'POSCAR file, plus the k.xyz extension.'
     case (1, 2)
        if ( i == 1 ) then
           write(arg(2), '(A)') trim(arg(1)) // '.k.xyz'
        end if
        call kspace2xyz (arg(1), arg(2))
     end select

  case ('check_nndist')
     select case (i)
     case default
        print *, 'Usage: check_nndist POSCAR [tol]'
        print *, ''
        print *, 'Check whether some atoms in the POSCAR file are closer '
        print *, 'to each other than a specified tolerance.'
        print *, 'POSCAR: The name of the input file.'
        print *, 'tol: the tolerance, in Ångströms. If omitted, &
             &defaults to 1.0 Å.'
     case (1, 2)
        if (i == 1) then
           call check_nn_POSCAR (arg(1))
        else if (i == 2) then
           read(arg(2), *) tol
           call check_nn_POSCAR (arg(1), tol)
        end if
     end select

  case ('scgenerator')
     select case (i)
     case default
        print *, 'Usage: scgenerator POSCAR POSCAR2 POSCAR.out'
        print *, ''
        print *, 'Generate a supercell.'
        print *, 'POSCAR: The name of the input file.'
        print *, 'POSCAR2: The optional name of the input file containing the &
             &supercell lattice vectors. If omitted, asks the user.'
        print *, 'POSCAR.out: The name of the output file. If omitted, &
             &defaults to the input file plus the .out extension'
        print *, 'The input file should contain the coordinates of the &
             &atoms to be replicated in the supercell using a space-filling &
             &algorithm. Note that you cannot omit only one of the arguments; &
             &the allowed forms are thus with 1 or 3 arguments.'
     case (1)
        write (arg(3), '(A)') trim (arg(1)) // '.out'
        call sc_generator_io (arg(1), outfile=arg(3))
     case (3)
!        if (i == 2) then
!           write(arg(3), '(A)') trim (arg(1)) // '.out'
!        end if
        call sc_generator_io (arg(1), arg(2), arg(3))
     end select

  case ('dumpcoords')
     select case (i)
     case default
        print *, 'Usage: vasputil dumpcoords POSCAR'
        print *, ''
        print *, 'Dump cartesian coordinates suitable for importing into '
        print *, 'e.g. octave or matlab.'
     case (1)
        call dumpcoords (arg(1))
     end select

  case ('importcoords')
     select case (i)
     case default
        print *, 'Usage: vasputil importcoords POSCAR < coords'
        print *, ''
        print *, 'Import cartesian coordinates from stdin.'
     case (1)
        call importcoords (arg(1))
     end select

  case ('test')
     call runtest ()

  case ('-V', '--version')
     call print_version ()
     
  case default
     call print_usage ()

  end select

contains

  
  !****f* vasputil/print_version
  ! PURPOSE
  ! Print the version number of the program to the screen.
  !****
  subroutine print_version ()
    print '(A, //, A, /, A, /, A)', ' vasputil release 3.3', &
         ' Copyright (C) 2004, 2005, 2006 Janne Blomqvist.', &
         ' This is free software; see the source for copying conditions.  &
         &There is NO', &
         ' warranty; not even for MERCHANTABILITY or FITNESS FOR &
         &A PARTICULAR PURPOSE.'
  end subroutine print_version


  !****f* vasputil/print_usage
  ! PURPOSE
  ! Print the usage information for vasputil.
  !****
  subroutine print_usage ()
    print *, 'Usage: vasputil [commandname|options] [command options] [files]'
    print *, 'vasputil consists of many utilities, all in the ', &
         'same binary that behaves differently depending on the ', &
         'name by which it is called.'
    print *, 'vasputil options:'
    print '(T4, A, T8, A, T32, A)', '-V,', '--version', &
         'output version information.'
    print *, ''
    print *, 'Available utilities are:'
    print *, 'poscar2xyz: Convert a POSCAR file to a xyz format file.'
    print *, 'plane2atom: Calculate the average distance from a &
         &plane to a set of atoms.'
    print *, 'plane2layer: Calculate the average distance from a &
         &plane to a layer.'
    print *, 'atomsdistance: Calculate the distance between two atoms.'
    print *, 'atomsmoved: Check which atoms have moved from one &
         &POSCAR to another.'
    print *, 'interpolate: Interpolate the coordinates between &
         &two POSCAR files.'
    print *, 'normalize: Normalize the coordinates in a POSCAR file.'
    print *, 'unnormalize: Unnormalize coordinates in a POSCAR file.'
    print *, 'removeatoms: Remove atoms from a POSCAR file.'
    print *, 'kspace2xyz: Convert a POSCAR file to a xyz file with &
         &reciprocal coordinates.'
    print *, 'check_nndist: Check nearest neighbor distances.'
    print *, 'scgenerator: Generate a supercell.'
    print *, 'lockatoms: Lock atom coordinates in a supercell.'
    print *, 'xyz2poscar: Convert from XYZ format to POSCAR.'
    print *, 'direct2cartesian: Convert POSCAR from direct to cartesian &
         &coordinates.'
    print *, 'dumpcoords: Dump cartesian coordinates to stdout.'
    print *, 'importcoords: Import cartesian coordinates from stdout.'
    print *, 'test: Print some diagnostic information.'
    print *, ''
    print *, 'Run the utilities without arguments to get usage &
         &instructions.'
  end subroutine print_usage


  !****f* vasputil/runtest
  ! PURPOSE
  ! Print some diagnostics information and run testsuite (whenever I
  ! get around to actually implementing the testsuite).
  !****
  subroutine runtest ()
    use kind_params
    print *, 'Real kinds supported by compiler:'
    print *, 'Single     Double   Quad      max(dp,qp)'
    print '(1X,I2,8X,I2,8X,I2,8X,I2)', sp, dp, qp_preferred, qp
  end subroutine runtest

end program vasputil
