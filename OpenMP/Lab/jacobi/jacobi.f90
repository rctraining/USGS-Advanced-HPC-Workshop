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

module jacobi

        use kinds,    only : r_dp
        use err,      only : err_msg
!        use hdf_rw,   only : hcreate, hwdata, hclose
!        use hdf5,     only : hid_t
        use walltime, only : get_walltime

        implicit none

        private
        public jacobi_solver

contains

subroutine jacobi_solver(phi, niter, ierr, filename)

        implicit none
        real(kind=r_dp), contiguous, intent(inout) :: phi(:,:,:)
        integer,                     intent(in)    :: niter
        integer,                     intent(out)   :: ierr
        character(len=*), optional,  intent(in)    :: filename

        integer                      :: it, k, j, i     ! loop indexers
        integer                      :: imax, jmax, kmax
        real(kind=r_dp), parameter   :: one_sixth = 1.0_r_dp / 6.0_r_dp
        real(kind=r_dp)              :: s, e
        real(kind=r_dp), allocatable :: tmp(:,:,:)
        !real(kind=r_dp), pointer     :: grid(:,:,:) => null()

        !integer(kind=hid_t) :: f_id                     ! file identifier
        !character(len=256)  :: outfile                  ! output file name
        !character(len=256)  :: gname                    ! iteration group

        ierr = 0
        imax = size(phi, 1)
        jmax = size(phi, 2)
        kmax = size(phi, 3)

        allocate(tmp, source=phi, stat=ierr)
        if (ierr /= 0 ) then
                call err_msg('Unable to allocate temporary domain')
        end if

        ! Write the array
!        if (present(filename)) then
!                outfile = trim(filename)
!        else
!                outfile = 'js.h5'
!        end if
!        call hcreate(outfile, f_id, ierr)
!        it = 0
!        write(gname, '(I8)') it
!        call hwdata(f_id, phi, gname, ierr)


        call get_walltime(s)

        do it = 1, niter
            do k = 2, kmax-1
                do j = 2, jmax-1
                    do i = 2, imax-1
                        tmp(i,j,k) =(  phi(i-1,j,k) + phi(i+1,j,k)  &
                                     + phi(i,j-1,k) + phi(i,j+1,k)  &
                                     + phi(i,j,k-1) + phi(i,j,k+1)) &
                                    * one_sixth
                    end do
                end do
            end do
            phi = tmp
!            write(gname, '(I8)') it
!            call hwdata(f_id, phi, gname, ierr)
        end do

        call get_walltime(e)

        print *, 'Time elapsed: ', e - s

!        call hclose(f_id, ierr)

        if (allocated(tmp)) then
                deallocate(tmp)
        end if

end subroutine jacobi_solver


end module jacobi
