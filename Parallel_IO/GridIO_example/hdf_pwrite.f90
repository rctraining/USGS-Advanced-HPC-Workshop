! Program to use MPI_Cart and Parallel HDF5
!
program hdf_pwrite

        use mpi
        use hdf5
        use kinds, only : r_dp

        implicit none

        ! Local array size with halo
        integer, parameter :: g_N   = 20
        integer, parameter :: ndims = 2
        integer, parameter :: halo  = 1

        integer :: argc         ! Number of command line arguments
        integer :: ierr         ! Error status
        integer :: id           ! My rank/ID
        integer :: np           ! Number of processors
        integer :: iunit        ! File descriptor
        integer :: i,j          ! Loop indexers
        integer :: n(ndims)     ! Local N for i and j directions
        integer :: total(ndims) ! Local total dimension size

        ! MPI IO/Lustre file striping
        integer :: lcount       ! Lustre count size
        integer :: lsize        ! Lustre stripe size
        character(len=1024) :: clcount, clsize ! Strings of LFS

        integer :: info                 ! MPI IO Info
        integer :: m_dims(ndims)        ! MPI cart dims
        integer :: coords(ndims)        ! Co-ords of procs in the grid
        logical :: is_periodic(ndims)   ! Periodic boundary conditions
        logical :: reorder              ! Reorder the MPI structure
        integer :: MPI_COMM_2D          ! New communicator

        character(len=1024) :: filename
        integer(kind=hid_t) :: p_id, f_id, x_id, d_id, c_id
        integer(kind=hid_t) :: memspace, filespace
        ! Chunk sizes
        integer(kind=hsize_t) :: c_size(ndims)
        ! Local hyper slab info
        integer(kind=hsize_t) :: d_size(ndims), s_size(ndims), h_size(ndims), &
                                 stride(ndims), block(ndims)
        ! Global hyper slab info
        integer(kind=hsize_t) :: g_size(ndims), g_start(ndims)

        ! Local data array
        real(kind=r_dp), allocatable :: grid(:,:)

        argc = 0
        ierr = 0
        m_dims = (/ 0, 0/)
        is_periodic = .true.     ! Periodic
        reorder     = .false.    ! Not allowed to reorder

        call mpi_init(ierr)

        ! Set up the MPI cartesian topology
        call mpi_comm_size(MPI_COMM_WORLD, np, ierr)
        call mpi_dims_create(np, ndims, m_dims, ierr)

        call mpi_cart_create(MPI_COMM_WORLD, ndims, m_dims, is_periodic, &
                             reorder, MPI_COMM_2D, ierr)
        call mpi_comm_rank(MPI_COMM_2D, id, ierr)
        call mpi_cart_coords(MPI_COMM_2D, id, ndims, coords, ierr)

        if (id .eq. 0) then
                if (mod(g_N,np) .ne. 0) then
                        write(0,*) 'Must use divisiable number of procs.'
                        call mpi_abort(MPI_COMM_WORLD, 1, ierr)
                endif

                ! get the filename
                argc = command_argument_count()
                if (argc .lt. 1 ) then
                        write(0, *) 'Must supply a filename'
                        call exit(1)
                endif
                call get_command_argument(1, filename)
        endif

        ! Broadcast the filename
        call mpi_bcast(filename, len(filename), MPI_CHAR, 0, &
                       MPI_COMM_WORLD, ierr)

        ! Init the HDF5 library
        call h5open_f(ierr)

        ! Set a stripe count of 4 and a stripe size of 4MB
        lcount = 4
        lsize  = 4 * 1024 * 1024
        write(clcount, '(I4)') lcount
        write(clsize, '(I8)') lsize

        call mpi_info_create(info, ierr)
        call mpi_info_set(info, "striping_factor", trim(clcount), ierr)
        call mpi_info_set(info, "striping_unit", trim(clsize), ierr)

        ! Set up the access properties
        call h5pcreate_f(H5P_FILE_ACCESS_F, p_id, ierr)
        call h5pset_fapl_mpio_f(p_id, MPI_COMM_2D, info, ierr)

        ! Open the file
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, f_id, ierr, &
                         access_prp = p_id)
        if (ierr .ne. 0) then
                write(0,*) 'Unable to open: ', trim(filename), ': ', ierr
                call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        endif

        ! Generate our local matrix
        do i = 1, ndims
                n(i) = g_N / m_dims(i)
                total(i) = n(i) + (2 * halo)
        end do

        allocate(grid(1 - halo : n(1) + halo,  &
                      1 - halo : n(2) + halo), &
                 stat=ierr)
        if (ierr .ne. 0) then
                write(0,*) 'Unable to allocate local data array: ', ierr
                call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        end if

        grid = -99.99
        ! Init the local data
        do j = 1, n(2)
                do i = 1, n(1)
                        grid(i,j) = id
                enddo
        enddo

        ! Create the local memory space and hyperslab
        do i = 1, ndims
                d_size(i) = total(i)
                s_size(i) = n(i)
                h_size(i) = halo
                stride(i) = 1
                block(i)  = 1
        enddo

        call h5screate_simple_f(ndims, d_size, memspace, ierr)
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                   h_size, s_size, ierr,       &
                                   stride, block)

        ! Create the global file space and hyperslab
        g_size  = g_N
        do i = 1, ndims
                g_start(i) = n(i) * coords(i)
        enddo

        call h5screate_simple_f(ndims, g_size, filespace, ierr)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
                                   g_start, s_size, ierr,       &
                                   stride, block)

        ! Create a data chunking property
        c_size = 2
        call h5pcreate_f(H5P_DATASET_CREATE_F, c_id, ierr)
        call h5pset_chunk_f(c_id, ndims, c_size, ierr)
        ! Create the dataset id
        call h5dcreate_f(f_id, "/data", H5T_IEEE_F64LE, filespace, d_id, &
                         ierr, dcpl_id=c_id)


        ! Create a data transfer property
        call h5pcreate_f(H5P_DATASET_XFER_F, x_id, ierr)
        call h5pset_dxpl_mpio_f(x_id, H5FD_MPIO_COLLECTIVE_F, ierr)
        ! Write the data
        call h5dwrite_f(d_id, H5T_IEEE_F64LE, grid, s_size, ierr,      &
                        file_space_id=filespace, mem_space_id=memspace, &
                        xfer_prp=x_id)

        if (allocated(grid)) then
                deallocate(grid)
        endif

        ! Close everything and exit
        call h5dclose_f(d_id, ierr)
        call h5sclose_f(filespace, ierr)
        call h5sclose_f(memspace, ierr)
        call h5pclose_f(c_id, ierr)
        call h5pclose_f(x_id, ierr)
        call h5pclose_f(p_id, ierr)
        call h5fclose_f(f_id, ierr)
        call h5close_f(ierr)

        call mpi_finalize(ierr)
end program hdf_pwrite
