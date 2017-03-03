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

! Jacobian Solver
program js

!        use hdf_rw, only : hinit
        use kinds,  only : r_dp
        use err,    only : err_abort
        use args,   only : args_parse
        use domain, only : domain_init
        use jacobi, only : jacobi_solver


        implicit none

        integer            :: ierr
        integer            :: niter
        integer            :: imax, jmax, kmax
        real(kind=r_dp)    :: north, south, east, west, top, bottom
        real(kind=r_dp), allocatable, target :: phi(:,:,:)
        character(len=256) :: filename

        ierr = 0
        filename = ''

        call args_parse(imax, jmax, kmax, niter,  &
                        north, south, east, west, &
                        top, bottom, ierr, output=filename)
        if (ierr /= 0) then
                call err_abort('Unable to parse input arguments', ierr)
        end if

         ! Create a new matrix and set the boundary conditions
        call domain_init(phi, imax, jmax, kmax, ierr, &
                         north, south, east, west, top, bottom)
        if (ierr /= 0) then
                call err_abort('Unable to initialize matrix', ierr)
        end if

        ! Initalize the HDF5 routines
!        call hinit(ierr)
!        if (ierr /= 0) then
!                call err_abort('Unable to initialize HDF5', ierr)
!        end if

        ! Solve the Jacobi
        if (filename /= '') then
                call jacobi_solver(phi, niter, ierr, filename)
        else
                call jacobi_solver(phi, niter, ierr)
        end if
        if (ierr /= 0) then
                call err_abort('Unable to run the jacobian solver', ierr)
        end if

        ! Deallocate the matrix
        if (allocated(phi)) then
                deallocate(phi)
        end if

end program js
