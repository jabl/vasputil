!****h* vasputil/supercell_generator
! NAME
!   supercell_generator
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
! This module defines a supercell generator.
!****
module supercell_generator
  use conf
  use geometry
  use lattice
  use supercell_core
  use supercell_utils
  use poscar_io

  implicit none

contains

  !****f* supercell_generator/generate_supercell
  ! PURPOSE
  !   Generate a supercell.
  ! INPUTS
  ! cell - The supercell which contains the primitive translation
  !        vectors and the basis set which is to be replicated into the
  !        supercell.
  ! scell - The supercell which is to be created. The lattice vectors of
  !         this cell should already be set up before calling this procedure.
  ! NOTES
  ! The algorithm works by using the lattice vectors from "cell" as the basis
  ! and using 
  ! those to fill the supercell. After that, the atoms in "cell" are
  ! replicated into "scell".
  ! BUGS
  ! Note that the variable names somewhat confusingly refer to
  ! atoms, but they are really mostly referring to the primitive
  ! translation vectors.
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
    integer :: natoms, i, j, natoms_created, hcol, hitot, hlookup, &
         natoms_tot, k
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
    if (debug) then
       print *, 'Estimated ', natoms, ' atoms in the new supercell.'
    end if

    ! Allocate the ready queue
    call init_queue (rlist, 2*natoms)

    ! This algorithm is a variation of a graph traversal algorithm, e.g. 
    ! "visit all the nodes in a graph". In this case "node" =
    ! primitive cell, i.e.
    ! primitive cells in the supercell
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

    if (debug) then
       if (hcol /= 0) then
          print *, 'A total of ', hlookup, &
               ' hash lookups were made, of which ', hcol, &
               ' resulted in collisions. The collisions were handled by a &
               &total of ', &
               hitot, 'rehashing operations.'
       end if
       print *, 'Explored ', j, ' atoms.'
    end if

    ! Now all the primitive cells have been created in the graph. 
    ! Next we walk
    ! through the graph again and fill the supercell coordinate array.
    if (debug) then
       print *, 'But really, we have ', natoms_created, ' basis cells in the supercell.'
    end if
    natoms_tot = natoms_created * size (cell%atoms)
    call init_cell (scell, natoms_tot)
    !      scell%atoms%symbol = cell%atoms(1)%symbol
    scell%cartesian = .true.
    scell%relative = .false.
    call direct2Cartesian (cell)
    call rel2Act (cell)
    ! Reinitialize the ready list
    call init_queue (rlist, natoms_created)
    avcur => NULL()
    call insert (rlist, av)
    j = 1 ! j is the 'first empty' index for scell%atomCoords.
    do k = 1, size(cell%atoms)
       scell%atomCoords(:, j) = cell%lattice%a * &
            matmul (cell%lattice%t, av%acoords + cell%atomCoords(:,k))
       scell%atoms(j)%symbol = cell%atoms(k)%symbol
       j = j + 1
    end do
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
                do k = 1, size (cell%atoms)
                   scell%atomCoords(:, j) = cell%lattice%a * &
                        matmul (cell%lattice%t, &
                        avnext%acoords + cell%atomCoords(:,k))
                   scell%atoms(j)%symbol = cell%atoms(k)%symbol
                   j = j + 1
                end do
                call insert (rlist, avnext)
             end if
          end if
       end do
    end do
    call sort_cell (scell)
    print *, 'The supercell contains ', j-1, ' atoms.'

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
      integer, parameter :: gr = 2654435761_i8b/2
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
    function is_inside_cell (scell, pcell, coords, tol)
      type(supercell), intent(in) :: scell
      type(supercell), intent(in) :: pcell
      integer, intent(in) :: coords(3)
      real(wp), intent(in), optional :: tol
      logical :: is_inside_cell, init = .false.
      real(wp) :: zert, onet, tole
      type(supercell), save :: cell ! Scratch cell.

      if (.not. init) then
         call init_cell (cell, 1)
         cell%lattice = pcell%lattice
         init = .true.
      end if
      if (present (tol)) then
         tole = tol
      else
         tole = 100 * sqrt (epsilon (1.0_wp))
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


  !****f* supercell_generator/sc_generator_io
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
    ! Make all coords cartesian and actual
    call direct2cartesian (cell)
    call direct2cartesian (scell)
    call rel2act (cell)
    call rel2act (scell)
    ! Fill out the new supercell
    call generate_supercell (cell, scell, error_stop)
    ! Write the new supercell
    call cartesian2Direct (scell)
    call act2rel (scell)
    call write_POSCAR (scell, outfile)
  end subroutine sc_generator_io


end module supercell_generator
