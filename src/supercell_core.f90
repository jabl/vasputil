!****h* vasputil/supercell_core
! NAME
!   supercell_core
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


    !****f* supercell_core/checkCells
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
         cell%atomCoords(dir, atoms(i)) = cell%atomCoords(dir, atoms(i)) - 1.0_wp
      end do
    end subroutine unnormalize


    !****f* supercell_core/atomdist
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
    

    !****f* supercell_core/planetopoint
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


    !****f* supercell_core/check_nndist
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


    !****f* supercell_core/generate_supercell
    ! PURPOSE
    !   Generate a supercell.
    ! INPUTS
    ! cell - The supercell which contains the basis set of the new supercell.
    ! scell - The supercell which is to be created. The lattice vectors of
    !         this cell should already be set up before calling this procedure.
    ! NOTES
    ! The algorithm works by using the lattice vectors from "cell" as the basis
    ! and using 
    ! those to fill the supercell.
    ! BUGS
    ! None. The code is perfect.
    !****
    subroutine generate_supercell (cell, scell, err_proc)
      type(supercell), intent(inout) :: cell
      type(supercell), intent(inout) :: scell
      interface
         subroutine err_proc (msg)
           character(len=*), intent(in) :: msg
         end subroutine err_proc
      end interface
      !real(wp) :: arel
      integer :: natoms, i, j, natoms_created, hcol, hitot, hlookup
      type atom_edge
         type(atom_vertex), pointer :: v => NULL()
      end type atom_edge
      type atom_vertex
         ! Whether we have visited all vertexes connected to this one.
         logical :: explored = .false.
         ! Atom coordinates.
         !real(wp) :: acoords(3) = (/ 0.0_wp, 0.0_wp, 0.0_wp /) 
         integer :: acoords(3) = (/0, 0, 0/)
         ! List of atoms connected to this one via a single basis vector.
         ! In 3D space we have 6 degrees of freedom
         type(atom_edge), dimension(6) :: alist
      end type atom_vertex
      ! Queue type, implemented using a circular array.
      type queue_t
         type(atom_edge), dimension(:), pointer :: elems => NULL()
         integer :: head = 1, tail = 1
      end type queue_t
      type(atom_vertex), target :: av
      type(atom_vertex), pointer :: avcur => NULL(), & ! Current atom vertex
           avnext => NULL() ! Next atom vertex
      ! The ready list
      type(queue_t) :: rlist

      if (is_lin_dependent (cell%lattice%t)) then
         call err_proc ("ERROR: Lattice vectors are not linearly independent!")
      end if

      call latt_init (scell%lattice)
      ! Calculate the ratio between the lattice constants
      !arel = scell%lattice%a / cell%lattice%a
      !print *, scell%lattice%omega, cell%lattice%omega
      ! Estimate for the number of atoms in the supercell
      natoms = ceiling (size (cell%atoms) * scell%lattice%omega / cell%lattice%omega)
      print *, 'Estimated ', natoms, ' atoms in the supercell.'

      ! Allocate the ready queue
      call init_queue (rlist, 2*natoms)

      ! This algorithm is a variation of a graph traversal algorithm, e.g. 
      ! "visit all the nodes in a graph". In this case "node" = atom, i.e. 
      ! atoms
      ! are created by starting from origo and adding the integer multiples of
      ! the basis vectors and checking that they are inside the bounding box
      ! created by the lattice vectors of scell.

      ! We do a breadth first search.
      call insert (rlist, av)
      natoms_created = 1
      hcol = 0
      hitot = 0
      hlookup = 0
      j = 1
      av%explored = .true.
      do
         call remove (rlist, avcur)
         if (.not. associated (avcur)) then
            exit
         end if
         call create_atom_children (avcur)
         do i = 1, size (avcur%alist)
            avnext => avcur%alist(i)%v
            if (associated (avnext)) then
               if (.not. avnext%explored) then
                  avnext%explored = .true.
                  j = j + 1
                  call insert (rlist, avnext)
               end if
            end if
         end do
      end do

      if (hcol /= 0) then
         print *, 'A total of ', hlookup, &
              ' hash lookups were made, of which ', hcol, &
              ' resulted in collisions. The collisions were handled by a &
              &total of ', &
              hitot, 'rehashing operations.'
      end if
      print *, 'Explored ', j, ' atoms.'

      ! Now all the atoms have been created in the graph. Next we walk
      ! through the graph again and fill the supercell coordinate array.
      print *, 'But really, we have ', natoms_created, ' atoms in the supercell.'
      call init_cell (scell, natoms_created)
      scell%atoms%symbol = cell%atoms(1)%symbol
      scell%cartesian = .true.
      scell%relative = .false.
      ! Reinitialize the ready list
      call init_queue (rlist, natoms_created)
      avcur => NULL()
      call insert (rlist, av)
      j = 1
      scell%atomCoords(:, j) = cell%lattice%a * matmul (cell%lattice%t, av%acoords)
      av%explored = .false.
      do
         call remove (rlist, avcur)
         if (.not. associated (avcur)) then
            exit
         end if
         do i = 1, size (avcur%alist)
            avnext => avcur%alist(i)%v
            if (associated (avnext)) then
               if (avnext%explored) then
                  avnext%explored = .false.
                  j = j + 1
                  scell%atomCoords(:, j) = cell%lattice%a * matmul (cell%lattice%t, avnext%acoords)
                  call insert (rlist, avnext)
               end if
            end if
         end do
      end do
      print *, 'Coordinates for ', j, ' atoms were copied to scell.'

      contains

        !****f* generate_supercell/hash
        ! PURPOSE
        !   The hash function used for atom coordinates.
        ! TODO
        !   Use a better hash function, this one creates lots of collisions.
        ! ARGUMENTS
        !   coords    - Input coordinates to hash.
        !   tsize     - Size of the hash table.
        !****
        pure function hash (coords, tsize)
          integer, intent(in) :: coords(3), tsize
          integer :: hash
          integer :: i, j, k, hkey
          ! Golden ratio of a signed 32-bit integer
          integer, parameter :: gr = 2654435761/2
          i = coords(1) * gr
          j = not (coords(2)) * gr
          k = ishftc (coords(3), 15) * gr
          ! \sum_{l}(i+j+k)^l  mod tsize
          !hash = modulo (i + j + k + (i + j + k)**2 &
          !     + (i + j + k)**3, tsize)
          ! Golden ratio of 2^32: 2654435761
          hkey = ieor (i, ior (j, k))
          hash = modulo (hkey, tsize) + 1
        end function hash


        !****f* generate_supercell/create_atom
        ! PURPOSE
        !   Given a new atom, initialize it by creating the children.
        !****
        subroutine create_atom_children (a)
          type(atom_vertex) :: a
          integer :: i, j, k, hval, coords(3), ind(3), hi
          ! Hash table for fast coordinate based lookup
          type(atom_edge), allocatable, save :: htable(:)
          
          if (.not. allocated (htable)) then
             allocate(htable(2*natoms + 1)) ! 2* is usually sparse enough.
             ! Insert the origin node
             hval = hash ((/0, 0, 0/), size (htable))
             htable(hval)%v => av
          end if

          k = 1
          do i = -1, 1, 2
             do j = 1, 3
                ind = 0
                ind(j) = i
                ! Check that we're inside the box
                coords = (/ a%acoords(1) + ind(1), a%acoords(2) + ind(2), &
                     a%acoords(3) + ind(3)/)
                if (.not. is_inside_cell (scell, cell, coords)) then
                   cycle
                end if
                hval = hash (coords, size (htable))
                !print *, i, j, hval, coords
                ! Check if the atom has already been created, in which case
                ! we just set the value of the pointer.
                hi = 1
                do
                   if (.not. associated (htable(hval)%v)) then
                      ! Create a new atom.
                      allocate(a%alist(k)%v)
                      natoms_created = natoms_created + 1
                      ! .. and update the hash table.
                      htable(hval)%v => a%alist(k)%v
                      ! Update the coords
                      a%alist(k)%v%acoords = coords
                      hlookup = hlookup + 1
                      exit
                   else if (htable(hval)%v%acoords(1) /= coords(1) .or. &
                        htable(hval)%v%acoords(2) /= coords(2) .or. &
                        htable(hval)%v%acoords(3) /= coords(3) ) then
                      ! Hash collision
                      if (hi == 1) then
                         hcol = hcol + 1
                      end if
                      hval = modulo (hval + 1 * hi**2, size (htable)) + 1
                      hi = hi + 1
                      !hlookup = hlookup + 1
                   else
                      ! We hit an existing atom that is not a hash collison,
                      ! so we point to it.
                      hitot = hitot + hi - 1
                      a%alist(k)%v => htable(hval)%v
                      hlookup = hlookup + 1
                      exit
                   end if
                end do
                k = k + 1
             end do
          end do
          !print *, 'Needed ', hitot, ' rehashing ops.'
        end subroutine create_atom_children


        !****f* generate_supercell/is_inside_cell
        ! PURPOSE
        ! Check that a coordinate point is inside the supercell boundary.
        ! NOTES
        ! Calculate the crossing number (Jordan curve theorem). Or rather,
        ! just check that direct coordinates are \in [0,1].
        !****
        function is_inside_cell (scell, cell, coords, tol)
          type(supercell), intent(in) :: scell
          type(supercell), intent(inout) :: cell
          integer, intent(in) :: coords(3)
          real(wp), intent(in), optional :: tol
          logical :: is_inside_cell
          real(wp) :: zert, onet, tole

          if (present (tol)) then
             tole = tol
          else
             tole = 1000 * sqrt (epsilon (1.0_wp))
          end if
          zert = 0.0_wp - tole
          onet = 1.0_wp - tole
          ! Get the coordinate in cartesian.
          cell%cartesian = .false.
          cell%relative = .true.
          cell%atomCoords(:,1) = (/real (coords(1), wp), &
               real (coords(2), wp), real (coords(3), wp) /)
          call direct2Cartesian (cell)
          call rel2Act (cell)
          ! Convert cartesian coordinate to scell direct coord.
          call cart2Direct (scell%lattice, cell%atomCoords)
          cell%atomCoords = cell%atomCoords / scell%lattice%a

          ! Finally, check whether coordinate is inside scell
          if (cell%atomCoords(1, 1) >= zert .and. &
               cell%atomCoords(1, 1) < onet .and. &
               cell%atomCoords(2, 1) >= zert .and. &
               cell%atomCoords(2, 1) < onet .and. &
               cell%atomCoords(3, 1) >= zert .and. &
               cell%atomCoords(3, 1) < onet) then
             is_inside_cell = .true.
          else
             is_inside_cell = .false.
          end if
        end function is_inside_cell


        !****f* generate_supercell/insert
        ! PURPOSE
        ! Insert an element into the ready queue.
        !****
        subroutine insert (queue, element)
          type(queue_t), intent(inout) :: queue
          type(atom_vertex), intent(in), target :: element

          queue%elems(queue%tail)%v => element
          queue%tail = modulo (queue%tail + 1, size (queue%elems))
        end subroutine insert

        
        !****f* generate_supercell/remove
        ! PURPOSE
        ! Remove an element from the queue (FIFO).
        !****
        subroutine remove (queue, element)
          type(queue_t), intent(inout) :: queue
          type(atom_vertex), pointer :: element

          element => queue%elems(queue%head)%v
          queue%elems(queue%head)%v => NULL()
          queue%head = modulo (queue%head + 1, size (queue%elems))
        end subroutine remove


        !****f* generate_supercell/init_queue
        ! PURPOSE
        ! Initialize a queue.
        !****
        subroutine init_queue (queue, qsize)
          type(queue_t), intent(inout) :: queue
          integer, intent(in) :: qsize
          integer :: i, j

          if (associated (queue%elems)) then
             j = 0
             do i = 1, size (queue%elems)
                if (associated (queue%elems(i)%v)) then
                   j = j + 1
                end if
             end do
             if (j /= 0) then
                print *, 'WARNING (init_queue): Associated elements in queue: ', j
             end if
             deallocate (queue%elems)
          end if
          allocate (queue%elems(qsize))
          queue%head = 1
          queue%tail = 1
        end subroutine init_queue


        !****f* generate_supercell/destroy_scgenerator
        ! PURPOSE
        ! Destructor procedure for the graph.
        ! BUGS
        ! So far it does (almost) nothing.
        ! NOTES
        ! Should do a depth first traversal of the graph and deallocate.
        !****
        subroutine destroy_scgenerator ()
          deallocate (rlist%elems)
        end subroutine destroy_scgenerator

    end subroutine generate_supercell


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

end module supercell_core
