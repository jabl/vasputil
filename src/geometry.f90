!****h* vasputil/geometry
! NAME
! geometry
! COPYRIGHT
!  Copyright (c) 2004, 2005, 2008 Janne Blomqvist

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
! Some geometric shapes and operations.
!****

module geometry
  use conf
  use f95_lapack

  implicit none

  !****t* geometry/plane3
  ! PURPOSE
  ! Plane in a 3d space.
  !****
  type plane3
     real(wp), dimension(3) :: normal
     real(wp) :: distOrigo
  end type plane3

  !****** geometry/.cross.
  ! PURPOSE
  ! Cross product operator.
  !****
  interface operator(.cross.)
     module procedure cross
  end interface

  !****** geometry/new
  ! PURPOSE
  ! Generic interface for creating a new plane.
  !****
  interface new
     module procedure new_pn, new_p3
  end interface

  contains

    !****f* geometry/new_pn
    ! PURPOSE
    ! Initialize a plane, using a point on the plane and the normal.
    !****
    subroutine new_pn(plane, point, normal)
      type(plane3), intent(out) :: plane
      real(wp), dimension(3), intent(in) :: point, normal
      plane%normal = normal
      plane%normal = plane%normal/norm(plane%normal)
      plane%distOrigo = dot_product(plane%normal, point)
    end subroutine new_pn


    !****f* geometry/new_p3
    ! PURPOSE
    ! Initialize a plane, using three points which lie in the plane.
    !****
    subroutine new_p3(plane, point1, point2, point3)
      type(plane3), intent(out) :: plane
      real(wp), dimension(3), intent(in) :: point1, point2, point3
      real(wp), dimension(3) :: v1, v2
      v1 = point2 - point1
      v2 = point3 - point2
      plane%normal = v1 .cross. v2
      plane%normal = plane%normal/norm(plane%normal)
      plane%distOrigo = dot_product(plane%normal, point2)
    end subroutine new_p3


    !****f* geometry/distanceFrom
    ! PURPOSE
    ! Distance between a plane and a point.
    !****
    pure function distanceFrom(plane, point)
      type(plane3), intent(in) :: plane
      real(wp), dimension(3), intent(in) :: point
      real(wp) :: distanceFrom
      
      distanceFrom = abs(dot_product(plane%normal, point) - plane%distOrigo)
    end function distanceFrom


    !****f* geometry/cross
    ! PURPOSE
    ! Cross product between two 3d vectors.
    !****
    pure function cross(a, b)
      real(wp), dimension(3), intent(in) :: a,b
      real(wp), dimension(3) :: cross

      cross = (/ a(2) * b(3) - a(3) * b(2), &
           a(3) * b(1) - a(1) * b(3), &
           a(1) * b(2) - a(2) * b(1) /)
    end function cross


    !****f* geometry/norm
    ! PURPOSE
    ! The length of a 3d vector.
    !****
    pure function norm(a)
      real(wp), dimension(3), intent(in) :: a
      real(wp) :: norm
      norm = sqrt(dot_product(a,a))
    end function norm


    !****f* geometry/pl_intersect
    ! PURPOSE
    !   Calcultes the point of intersection between a line and a plane
    ! NOTES
    ! See e.g. http://astronomy.swin.edu.au/~pbourke/geometry/planeline/
    !****
    subroutine pl_intersect (plane, linep1, linep2, intersect, parallel)
      type(plane3), intent(in) :: plane
      real(wp), dimension(3), intent(in) :: linep1, linep2
      real(wp), optional, intent(out) :: intersect(3)
      logical, optional, intent(out) :: parallel
      real(wp) :: u, uden, unom
      ! If the denominator is 0, the line is parallel to the plane.
      uden = dot_product (plane%normal, linep2 - linep1)
      if (present (parallel)) then
         if (abs (uden) < epsilon (uden)) then
            parallel = .true.
         else
            parallel = .false.
         end if
      end if
      unom = dot_product (plane%normal, plane%distOrigo * plane%normal - linep1)
      u = unom / uden
      ! The point of intersection
      if (present (intersect)) then
         intersect = linep1 + u * (linep2 - linep1)
      end if
    end subroutine pl_intersect


    !****f* geometry/is_lin_dependent
    ! PURPOSE
    ! Check if a set of vectors are linearly dependent
    !****
    function is_lin_dependent (a)
      real(wp), dimension(:,:), intent(in) :: a
      logical :: is_lin_dependent
      real(wp) :: s(3), icond, a2(size (a, 1), size (a, 2))

      a2 = a
      ! Perform SVD via LAPACK95
      call la_gesdd (a2, s)

      ! Condition number of the matrix is the ratio between 
      ! the largest to the smallest singular value. If this
      ! approaches infinity, the matrix is ill-conditioned, 
      ! i.e. the input vectors were not linearly independent.
      ! For numerical reasons calculate the inverse of the 
      ! condition number.
      icond = minval (s) / maxval (s)
      if (icond < epsilon (icond)) then
         is_lin_dependent = .true.
      else
         is_lin_dependent = .false.
      end if
    end function is_lin_dependent

end module geometry
