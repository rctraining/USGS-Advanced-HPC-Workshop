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
!
! Module to read and write a dataset to a HDF5
!
module hdf_rw

        use hdf5,  only : hid_t, hsize_t,                               &
                          H5P_DEFAULT_F,                                &
                          H5F_ACC_TRUNC_F, H5F_ACC_RDWR_F,              &
                          H5T_NATIVE_DOUBLE,                            &
                          h5open_f, h5fcreate_f, h5fopen_f, h5fclose_f, &
                          h5gcreate_f, h5gclose_f,                      &
                          h5screate_simple_f, h5sclose_f,               &
                          h5dcreate_f, h5dwrite_f, h5dclose_f
        use kinds, only : r_dp
        use err,   only : err_msg

        implicit none

        private
        public  :: hinit,   &
                   hcreate, &
                   hopen,   &
                   hrdata,  &
                   hwdata,  &
                   hclose

contains

subroutine hinit(ierr)

        implicit none
        integer,             intent(out) :: ierr      ! error status

        call h5open_f(ierr)
        if (ierr /= 0) then
                call err_msg('Unable to initalise HDF5')
        end if

end subroutine hinit

subroutine hcreate(filename, f_id, ierr)

        implicit none

        character(len=*),    intent(in)  :: filename  ! file name
        integer(kind=hid_t), intent(out) :: f_id      ! file id
        integer,             intent(out) :: ierr      ! error status

        ! Create the file
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f_id, ierr)
        if (ierr /= 0) then
                call err_msg('Unable to create a HDF5 file: ' &
                         // trim(filename))
                return
        end if

end subroutine hcreate

subroutine hopen(filename, f_id, ierr)

        implicit none

        character(len=*),    intent(in)  :: filename  ! file name
        integer(kind=hid_t), intent(out) :: f_id      ! file id
        integer,             intent(out) :: ierr      ! error status

        ! Open the file
        call h5fopen_f(filename, H5F_ACC_RDWR_F, f_id, ierr)
        if (ierr /= 0) then
                call err_msg('Unable to create a HDF5 file: ' &
                         // trim(filename))
                return
        end if

end subroutine hopen

subroutine hclose(f_id, ierr)

        implicit none

        integer(kind=hid_t), intent(in)  :: f_id
        integer,             intent(out) :: ierr

        call h5fclose_f(f_id, ierr)
        if (ierr /= 0) then
                call err_msg('Unable to close HDF5 file')
                return
        end if

end subroutine hclose

subroutine hwdata(f_id, data, group, ierr)
        implicit none

        integer(kind=hid_t), intent(in)  :: f_id
        real(kind=r_dp),     intent(in)  :: data(:,:,:)
        character(len=*),    intent(in)  :: group
        integer,             intent(out) :: ierr

        character(len=256), parameter    :: dname = 'data'
        integer(kind=hid_t)              :: g_id      ! group identifier
        integer(kind=hid_t)              :: dset_id   ! dataset identifier
        integer(kind=hid_t)              :: dspace_id ! dataspace identifier

        integer, parameter               :: rank = 3
        integer(kind=hsize_t), dimension(3) :: dims

        ierr = 0
        dims(1) = size(data, 1)
        dims(2) = size(data, 2)
        dims(3) = size(data, 3)

        call h5gcreate_f(f_id, group, g_id, ierr)

        call h5screate_simple_f(rank, dims, dspace_id, ierr)

        call h5dcreate_f(g_id, dname, H5T_NATIVE_DOUBLE, dspace_id, &
                       dset_id, ierr, H5P_DEFAULT_F, H5P_DEFAULT_F, &
                       H5P_DEFAULT_F)

        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, ierr)

        call h5dclose_f(dset_id, ierr)

        call h5sclose_f(dspace_id, ierr)



end subroutine hwdata

subroutine hrdata(f_id, data, ierr)
        implicit none

        integer(kind=hid_t), intent(in)  :: f_id
        integer,             intent(in)  :: data
        integer,             intent(out) :: ierr

        ierr = 0
end subroutine hrdata

!subroutine write_attr()
!end subroutine write_attr
!
!        ! Local 4x4 with a 1x1 halo
!        integer, parameter :: ndims = 2
!        integer, parameter :: N     = 4000
!        integer, parameter :: halo  = 1
!
!        integer :: argc
!        integer :: ierr         ! Error status
!        integer :: id           ! My rank/ID
!        integer :: np           ! Number of processors
!        integer :: iunit        ! File descriptor
!        integer :: i,j          ! Loop indexers
!        integer :: total        ! Total dimension size
!        integer :: lcount       ! Luster count size
!        integer :: lsize        ! Lustre stripe size
!        character(len=1024) :: clcount, clsize ! Strings of LFS
!
!        integer :: info                 ! MPI IO Info
!        integer :: m_dims(ndims)        ! MPI cart dims
!        integer :: coords(ndims)        ! Co-ords of procs in the grid
!        logical :: is_periodic(ndims)   ! Periodic boundary conditions
!        logical :: reorder              ! Reorder the MPI structure
!        integer :: MPI_COMM_2D          ! New communicator
!
!        integer(KIND=MPI_OFFSET_KIND) :: offset
!
!        character(len=1024) :: filename
!        integer(kind=hid_t) :: p_id, f_id, x_id, d_id
!        integer(kind=hid_t) :: memspace, filespace
!        ! Local hyper slab info
!        integer(kind=hsize_t) :: d_size(ndims), s_size(ndims), h_size(ndims), &
!                                 stride(ndims), block(ndims)
!        ! Global hyper slab info
!        integer(kind=hsize_t) :: g_size(ndims), g_start(ndims)
!
!        real(kind=r_dp), allocatable :: ld(:,:)
!        ! Timing vars
!        real(kind=r_dp) :: s, e, dt, mdt
!
!        argc = 0
!        ierr = 0
!        offset = 0
!        m_dims = (/ 0, 0/)
!        is_periodic = .false.      ! Non-periodic
!        reorder     = .false.      ! Not allowed to reorder
!
!        call mpi_init(ierr)
!
!        ! Set up the MPI cartesian topology
!        call mpi_comm_size(MPI_COMM_WORLD, np, ierr)
!        call mpi_dims_create(np, ndims, m_dims, ierr)
!
!        call mpi_cart_create(MPI_COMM_WORLD, ndims, m_dims, is_periodic, &
!                             reorder, MPI_COMM_2D, ierr)
!        call mpi_comm_rank(MPI_COMM_2D, id, ierr)
!        call mpi_cart_coords(MPI_COMM_2D, id, ndims, coords, ierr)
!
!        if (id .eq. 0) then
!                if (mod(N,np) .ne. 0) then
!                        write(0,*) 'Must use divisiable number of procs.'
!                        call mpi_abort(MPI_COMM_WORLD, 1, ierr)
!                end if
!
!                ! get the filename
!                argc = iargc()
!                if (argc .lt. 1 ) then
!                        write(0, *) 'Must supply a filename'
!                        call exit(1)
!                end if
!                call get_command_argument(1, filename)
!        end if
!
!        ! Broadcast the filename
!        call mpi_bcast(filename, len(filename), MPI_CHAR, 0, &
!                       MPI_COMM_WORLD, ierr)
!
!        ! Init the HDF5 library
!        call h5open_f(ierr)
!
!        ! Set a stripe count of 4 and a stripe size of 4MB
!        lcount = 4
!        lsize  = 4 * 1024 * 1024
!        write(clcount, '(I4)') lcount
!        write(clsize, '(I8)') lsize
!
!        call mpi_info_create(info, ierr)
!        call mpi_info_set(info, "striping_factor", trim(clcount), ierr)
!        call mpi_info_set(info, "striping_unit", trim(clsize), ierr)
!
!        ! Set up the access properties
!        call h5pcreate_f(H5P_FILE_ACCESS_F, p_id, ierr)
!        call h5pset_fapl_mpio_f(p_id, MPI_COMM_2D, info, ierr)
!
!        ! Open the file
!        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f_id, ierr, &
!                         access_prp = p_id)
!        if (ierr .ne. 0) then
!                write(0,*) 'Unable to open: ', trim(filename), ': ', ierr
!                call mpi_abort(MPI_COMM_WORLD, 1, ierr)
!        end if
!
!        ! Generate our 4x4 matrix with a 1x1 halo
!        total = N + 2 * halo
!        allocate(ld(0:total-1, 0:total-1))
!
!        ld = -99.99
!        ! init the local data
!        do j = 1, N
!                do i = 1, N
!                        ld(i,j) = (i - 1 + (j-1)*N)
!                enddo
!        enddo
!
!        ! Create the local memory space and hyperslab
!        do i = 1, ndims
!                d_size(i) = total
!                s_size(i) = N
!                h_size(i) = halo
!                stride(i) = 1
!                block(i)  = 1
!        enddo
!
!        call h5screate_simple_f(ndims, d_size, memspace, ierr)
!        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
!                                   h_size, s_size, ierr,       &
!                                   stride, block)
!
!        ! Create the global file space and hyperslab
!        do i = 1, ndims
!                g_size(i)  = N * m_dims(i)
!                g_start(i) = N * coords(i)
!        enddo
!
!        call h5screate_simple_f(ndims, g_size, filespace, ierr)
!        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
!                                   g_start, s_size, ierr,       &
!                                   stride, block)
!
!        ! Create a data transfer property
!        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
!        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
!
!        ! Create the dataset id
!        call h5dcreate_f(f_id, "/data", H5T_IEEE_F64LE, filespace, d_id, &
!                         ierr)
!
!
!        ! Write the data
!        call get_walltime(s)
!        call h5dwrite_f(d_id, H5T_NATIVE_DOUBLE, ld, s_size, ierr,      &
!                        file_space_id=filespace, mem_space_id=memspace, &
!                        xfer_prp=x_id)
!        call get_walltime(e)
!
!        dt = e - s
!        call mpi_reduce(dt, mdt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_2D, ierr)
!
!        if (id .eq. 0) then
!                write(6,*) mdt / np
!        end if
!
!        if (allocated(ld)) then
!                deallocate(ld)
!        end if
!
!        ! Close everything and exit
!        call h5dclose_f(d_id, ierr)
!        call h5sclose_f(filespace, ierr)
!        call h5sclose_f(memspace, ierr)
!        call h5pclose_f(x_id, ierr)
!        call h5pclose_f(p_id, ierr)
!        call h5fclose_f(f_id, ierr)
!        call h5close_f(ierr)
!
!        call mpi_finalize(ierr)
end module hdf_rw
