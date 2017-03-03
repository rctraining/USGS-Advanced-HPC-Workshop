!
! Copyright (C) 2016  Timothy Brown
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!

!
! Module for domain decomposition.
!
! In this module the following definitions are used for the grid
! definition, sides and neighbours.
!
!
!
!                      above
!                        ^    up                          k
!                        |    /                           ^  j
!                        |   /                            |  ^
!               +--------|---------+                      | /
!              /|        |        /|                      |/
!             / |                / |                      +----> i
!            +------------------+  |
!    left    |  |               |  |     right
!    <-------|- |               | ---------->
!            |  |      /        |  |
!            |  +-----/---------|--+
!            | /     /          | /
!            |/     /   |       |/
!            +-----/------------+
!                 /     |
!                /      |
!              down     v
!                      below
!
!
!
!


module domain

        use kinds, only : r_dp
        use err,   only : err_msg

        implicit none

        private

        public  :: domain_init, domain_bc

contains

! Initalize the domain/grid
subroutine domain_init(grid, imax, jmax, kmax, ierr, &
                       north, south, east, west, top, bottom)

        implicit none

        real(kind=r_dp), allocatable, intent(inout) :: grid(:,:,:)
        integer,                      intent(in)    :: imax
        integer,                      intent(in)    :: jmax
        integer,                      intent(in)    :: kmax
        integer,                      intent(out)   :: ierr
        real(kind=r_dp), optional,    intent(in)    :: north
        real(kind=r_dp), optional,    intent(in)    :: south
        real(kind=r_dp), optional,    intent(in)    :: east
        real(kind=r_dp), optional,    intent(in)    :: west
        real(kind=r_dp), optional,    intent(in)    :: top
        real(kind=r_dp), optional,    intent(in)    :: bottom

        ierr = 0

        if (allocated(grid)) then
                call err_msg('Domain is already initialized')
                ierr = 1
                return
        end if

        allocate(grid(imax, jmax, kmax), stat=ierr)
        if (ierr /= 0) then
                call err_msg('Unable to allocate new grid')
                return
        end if
        grid = 0.0_r_dp

        call domain_bc(grid, ierr, north, south, east, west, top, bottom)
        if (ierr /= 0) then
                call err_msg('Unable to set boundary conditions on grid')
                return
        end if

end subroutine domain_init

! Set the initial boundary conditions
subroutine domain_bc(grid, ierr, north, south, east, west, top, bottom)

        implicit none

        real(kind=r_dp),           intent(inout) :: grid(:,:,:)
        integer,                   intent(out)   :: ierr
        real(kind=r_dp), optional, intent(in)    :: north
        real(kind=r_dp), optional, intent(in)    :: south
        real(kind=r_dp), optional, intent(in)    :: east
        real(kind=r_dp), optional, intent(in)    :: west
        real(kind=r_dp), optional, intent(in)    :: top
        real(kind=r_dp), optional, intent(in)    :: bottom

        integer :: imax
        integer :: jmax
        integer :: kmax
        real(kind=r_dp) :: r_north
        real(kind=r_dp) :: r_south
        real(kind=r_dp) :: r_east
        real(kind=r_dp) :: r_west
        real(kind=r_dp) :: r_top
        real(kind=r_dp) :: r_bottom

        ierr = 0
        imax = size(grid, 1)
        jmax = size(grid, 2)
        kmax = size(grid, 3)

        if (present(west)) then
                r_west = west
        else
                r_west = 0.0_r_dp
        end if

        if (present(east)) then
                r_east = east
        else
                r_east = 0.0_r_dp
        end if

        if (present(south)) then
                r_south = south
        else
                r_south = 0.0_r_dp
        end if
        if (present(north)) then
                r_north = north
        else
                r_north = 0.0_r_dp
        end if
        if (present(bottom)) then
                r_bottom = bottom
        else
                r_bottom = 0.0_r_dp
        end if
        if (present(top)) then
                r_top = top
        else
                r_top = 0.0_r_dp
        end if

        ! Faces
        grid(1,        2:jmax-1, 2:kmax-1) = r_west
        grid(imax,     2:jmax-1, 2:kmax-1) = r_east
        grid(2:imax-1, 2:jmax-1, 1       ) = r_south
        grid(2:imax-1, 2:jmax-1, kmax    ) = r_north
        grid(2:imax-1, 1,        2:kmax-1) = r_bottom
        grid(2:imax-1, jmax,     2:kmax-1) = r_top

        ! Edge swaths
        grid(2:imax-1, 1,        1       ) = max(r_bottom, r_south)
        grid(1,        2:jmax-1, 1       ) = max(r_west, r_south)
        grid(imax,     2:jmax-1, 1       ) = max(r_east, r_south)
        grid(2:imax-1, jmax,     1       ) = max(r_south, r_top)
        grid(2:imax-1, 1,        kmax    ) = max(r_bottom, r_north)
        grid(1,        2:jmax-1, kmax    ) = max(r_west, r_north)
        grid(imax,     2:jmax-1, kmax    ) = max(r_north, r_east)
        grid(2:imax-1, jmax,     kmax    ) = max(r_north, r_top)
        grid(1,        1,        2:kmax-1) = max(r_west, r_bottom)
        grid(imax,     1,        2:kmax-1) = max(r_east, r_bottom)
        grid(imax,     jmax,     2:kmax-1) = max(r_east, r_top)
        grid(1,        jmax,     2:kmax-1) = max(r_west, r_top)

        ! Corners
        grid(1,    1,    1   ) = max(r_west, max(r_south, r_bottom))
        grid(imax, 1,    1   ) = max(r_east, max(r_south, r_bottom))
        grid(1,    1,    kmax) = max(r_west, max(r_north, r_bottom))
        grid(imax, 1,    kmax) = max(r_east, max(r_north, r_bottom))
        grid(1,    jmax, 1   ) = max(r_west, max(r_south, r_top))
        grid(1,    jmax, kmax) = max(r_west, max(r_north, r_top))
        grid(imax, jmax, 1   ) = max(r_east, max(r_south, r_top))
        grid(imax, jmax, kmax) = max(r_east, max(r_top,   r_north))

end subroutine domain_bc

end module domain
