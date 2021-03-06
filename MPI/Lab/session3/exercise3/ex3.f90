!//////////////////////////////////////////////////////////
! Exercise : 2-D Diffusion Program
!            Parallelize this program by filling in the body of GHOST_ZONE_COMM
!            Based the parallelization on the 2-D domain decomposition by
!            communicating using the row and column communicators and their associated
!            row and column ranks.
!            NOTE:  It will probably help you to copy and modify 
!                   one of the ghost_zone_comm routines from session 1
PROGRAM MAIN
    USE MPI
    USE TIMING
    USE ISO_FORTRAN_ENV, ONLY : output_unit
    IMPLICIT NONE
    ! We use some objects of the Timer class defined in Timing.f90
    ! to profile our code.
    TYPE(Timer) :: calc_time, comm_time, init_time, wall_time
    REAL*8 :: calc_ratio, comm_ratio, init_ratio
    REAL*8, ALLOCATABLE :: var(:,:,:), wrk(:,:,:)
    INTEGER ::  nx,ny,nz,nt
    INTEGER :: i,j,k,n, ierr
    CHARACTER*120 :: msg
    !/////////////////////////////////////////////////////////////////////////
    !  MPI Process-Grid Variables
    INTEGER :: my_rank, ncpu ! global rank and number of MPI processes

    INTEGER :: my_row_rank, my_col_rank ! The rank within a row and column for each process respectively

    INTEGER :: nprow, npcol ! Number of process rows and columns respectively
                               ! 0 <= my_row_rank <= npcol
                               ! 0 <= my_col_rank <= nprow

    INTEGER :: col_comm, row_comm ! subcommunicators for rows and columns
                                   ! row_comm connects all processes in the same row
                                   ! col_comm connects all processes in the same column

    !/////////////////////////////////////////
    ! Variables related to load-balancing in x-direction.
    INTEGER :: my_imin, my_imax ! the range of x-indices that this rank owns
    INTEGER :: my_numx

    !/////////////////////////////////////////
    ! Variables related to load-balancing in y-direction.
    INTEGER :: my_jmin, my_jmax ! the range of y-indices that this column owns
    INTEGER :: my_numy

    !////////////////////////////////////////////////////
    ! Set some default values for the problem parameters.
    nx = 256   ! x-dimension of our 3-D array (global)
    ny = 256   ! y-dimension of our 3-D array (global)
    nz = 1     ! z-dimension of our 3-D array
    nt = 100   ! number of time steps to integrate over
    nprow = 2  ! number of process rows/number of ranks within a column
    npcol = 2  ! number of process columns/number of ranks within a row

    !////////////////////////////////////////////////////
    ! Check to see if the user has overridden
    ! these defaults at the command line.
    ! Command line calling syntax:
    ! ./ex5.gpu -nx 256 -ny 256 -nz 256 -nt 100 -nprow 2 -npcol 2

    !/////////////////////// Initialization
    CALL GRAB_ARGS()
    nz = 1  !hard-coded as 2-D for now.  Will expand to 3-D later.
    CALL INIT_COMM()

    CALL calc_time%init()
    CALL comm_time%init()
    CALL init_time%init()
    CALL wall_time%init()

    CALL wall_time%startclock()
    CALL init_time%startclock()

    CALL LOAD_BALANCE()

    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,'(a)')   ' ///////////////////////////'
        WRITE(output_unit,'(a)')   '  3-D Diffusion Equation   '
        WRITE(output_unit,'(a)')   '  --Problem Parameters-- '
        WRITE(output_unit,'(a,i0)')'    nx = ', nx
        WRITE(output_unit,'(a,i0)')'    ny = ', ny
        WRITE(output_unit,'(a,i0)')'    nz = ', nz
        WRITE(output_unit,'(a,i0)')'    nt = ', nt
        WRITE(output_unit,'(a,i0)')'  ncpu = ', ncpu
        WRITE(output_unit,'(a,i0)')' nprow = ', nprow
        WRITE(output_unit,'(a,i0)')' npcol = ', npcol
        WRITE(output_unit,'(a)')   ' ///////////////////////////'
        WRITE(output_unit,'(a)')   ' '
        WRITE(output_unit,'(a)')   ' Proccess Grid Layout: '
    ENDIF

    WRITE(msg,'(a,i0,a,i0,a,i0)') &
       ' row = ', my_col_rank, ' column = ', my_row_rank, ' rank = ', my_rank
    CALL ORDERED_PRINT(msg)


    !/////////////////////////////////////////////////////
    ! Initialize var and the work array
    ALLOCATE(var(my_imin-1:my_imax+1,my_jmin-1:my_jmax+1,1:nz))
    ALLOCATE(wrk(1:nx,1:ny,1:nz))
    
    var(:,:,:) = 0.0d0
    wrk(:,:,:)   = 0.0d0
    CALL INIT_ARR(var, 1.0d0, 2, 2, 2) ! product of single sine waves in each dimension


    CALL init_time%increment()


    !/////////////////////////////////////////////////////
    ! Evolve the system

    DO n = 1, nt
        IF (my_rank .eq. 0) THEN
            IF (MOD(n,100) .eq. 0) Write(output_unit,'(a,i5)')' Timestep: ',n
        ENDIF

        !* EDIT GHOST_ZONE_COMM to communicate ghost zones
        !  to left and right neighbor ranks.

        ! Communication
        CALL comm_time%startclock()
        IF (ncpu .gt. 1) CALL GHOST_ZONE_COMM()
        CALL comm_time%increment()


        ! Calculation
        CALL calc_time%startclock()
        IF (nz .eq. 1) THEN
            CALL Laplacian2D()
        ELSE
            CALL Laplacian3D()
        ENDIF
        CALL calc_time%increment()

    ENDDO

    CALL wall_time%increment

    init_ratio = init_time%elapsed/wall_time%elapsed*100
    calc_ratio = calc_time%elapsed/wall_time%elapsed*100
    comm_ratio = comm_time%elapsed/wall_time%elapsed*100

    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,'(a)')' '
        WRITE(output_unit,'(a)')' Complete!'
        WRITE(output_unit,'(a)')' '
        WRITE(output_unit,'(a)')' Measured Timings (from rank 0):'
        WRITE(output_unit,'(a)')' ...............................'
        WRITE(output_unit,'( a, ES14.4, a)') &
            '         Elapsed time: ', wall_time%elapsed, ' seconds.'
        WRITE(output_unit, fmt= '( a, ES14.4, a, F6.3,a)') &
            '  Initialization time: ', init_time%elapsed, ' seconds (',init_ratio,' %).'
        WRITE(output_unit, fmt= '( a, ES14.4, a, F6.3,a)') &
            '     Calculation time: ', calc_time%elapsed, ' seconds (',calc_ratio,' %).'
        WRITE(output_unit, fmt= '( a, ES14.4, a, F6.3,a)') &
            '   Communication time: ', comm_time%elapsed, ' seconds (',comm_ratio,' %).'
    ENDIF

    CALL MPI_FINALIZE(ierr)
    
CONTAINS

    SUBROUTINE GHOST_ZONE_COMM()
        IMPLICIT NONE

    END SUBROUTINE GHOST_ZONE_COMM

    SUBROUTINE ORDERED_PRINT(message)
        IMPLICIT NONE
        INTEGER :: i,j, ierr
        CHARACTER(*), INTENT(IN) :: message
        DO j = 0, nprow-1    ! Loop over each row
            DO i = 0, npcol-1 ! Loop over each column
                IF ((my_row_rank .eq. i) .and. (my_col_rank .eq. j) ) THEN
                    WRITE(output_unit,*) TRIM(msg)
                ENDIF
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            ENDDO
        ENDDO
    END SUBROUTINE ORDERED_PRINT

    SUBROUTINE INIT_COMM()
        IMPLICIT NONE
        INTEGER :: ierr
        !Initialize MPI  (ierr is an error flag with nonzero value if there's a problem)
        Call MPI_INIT( ierr )
        IF (ierr .ne. 0) STOP 'Error intializing MPI'

        !Find number of MPI processes executing this program (ncpu)
        Call MPI_Comm_size(MPI_COMM_WORLD, ncpu,ierr)
        IF (ierr .ne. 0) STOP 'Error finding ncpu'
        IF (ncpu .ne. (nprow*npcol)) THEN
            STOP 'Error!  ncpu MUST EQUAL nprow*npcol'
        ENDIF


        !Find this process's rank (my_rank; unique numeric identifier for each MPI process)
        Call MPI_Comm_rank(MPI_COMM_WORLD, my_rank,ierr)
        IF (ierr .ne. 0) STOP 'Error finding my_rank'

        Call ROW_COLUMN_INIT()    

        !WRITE(output_unit,'(a,i0,a,i0,a)')'Rank ',my_rank,' of ', ncpu, ' is active.'
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 

    END SUBROUTINE INIT_COMM


    SUBROUTINE ROW_COLUMN_INIT()
        IMPLICIT NONE
        INTEGER :: test, ierr
        my_row_rank = MOD(my_rank,npcol)
        my_col_rank = my_rank/npcol

        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, my_col_rank, my_rank, row_comm, ierr)

        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, my_row_rank, my_rank, col_comm, ierr)
    END SUBROUTINE ROW_COLUMN_INIT

    SUBROUTINE grab_args()
            IMPLICIT NONE
            INTEGER :: n                    ! Number of command-line arguments
            INTEGER :: i                    
            CHARACTER(len=1024) :: argname  ! Argument key
            CHARACTER(len=1024) :: val      ! Argument value

            n = command_argument_count()
            DO i=1,n,2
                    CALL get_command_argument(i, argname)
                    CALL get_command_argument(i+1, val)
                    SELECT CASE(argname)
                            CASE('-nprow')
                                    read(val, '(I8)') nprow
                            CASE('-npcol')
                                    read(val, '(I8)') npcol
                            CASE('-nx')
                                    read(val, '(I8)') nx
                            CASE('-ny')
                                    read(val, '(I8)') ny
                            CASE('-nz')
                                    read(val, '(I8)') nz
                            CASE('-nt')
                                    read(val, '(I8)') nt
                            CASE DEFAULT
                                    WRITE(output_unit,'(a)') ' '
                                    WRITE(output_unit,'(a)') &
                                    ' Unrecognized option: '// trim(argname)
                    END SELECT
            ENDDO


    END SUBROUTINE grab_args

    SUBROUTINE Laplacian2D()
        IMPLICIT NONE
        Real*8 :: one_fourth = 0.25d0
        INTEGER :: dims(3), ni, nj, nk
        INTEGER :: imin,imax, jmin,jmax

        imin = max(my_imin,2)    ! Avoid the boundaries, which are held fixed at zero
        imax = min(my_imax,nx-1)

        jmin = max(my_jmin,2)    ! Avoid the boundaries, which are held fixed at zero
        jmax = min(my_jmax,ny-1)


        DO j = jmin,jmax
            DO i = imin,imax
                wrk(i,j,1) =  &
                        var(i-1,j,1) + var(i+1,j,1) + &
                        var(i,j-1,1) + var(i,j+1,1)
            ENDDO
        ENDDO

        DO j = jmin,jmax
            DO i = imin,imax
                var(i,j,1) = wrk(i,j,1)*one_fourth
            ENDDO
        ENDDO

    END SUBROUTINE Laplacian2D

    SUBROUTINE Laplacian3D()
        IMPLICIT NONE
        Real*8 :: one_sixth
        INTEGER :: imin,imax,jmin,jmax

        one_sixth = 1.0d0/6.0d0

        imin = max(my_imin,2)    ! Avoid the boundaries, which are held fixed at zero
        imax = min(my_imax,nx-1)

        jmin = max(my_jmin,2)    ! Avoid the boundaries, which are held fixed at zero
        jmax = min(my_jmax,ny-1)

        DO k = 2, nz-1
            DO j = jmin,jmax
                DO i = imin,imax
                    wrk(i,j,k) =  &
                            var(i-1,j,k) + var(i+1,j,k) + &
                            var(i,j-1,k) + var(i,j+1,k) + & 
                            var(i,j,k-1) + var(i,j,k+1)
                ENDDO
            ENDDO
        ENDDO

        DO k = 2, nz-1
            DO j = jmin,jmax
                DO i = imin,imax
                    var(i,j,k) = wrk(i,j,k)*one_sixth
                ENDDO
            ENDDO
        ENDDO


    END SUBROUTINE Laplacian3D


    SUBROUTINE LOAD_BALANCE()
        IMPLICIT NONE
        INTEGER :: modcheck
        !FIND (imin,imax) & (jmin,jmax) for this process

        ! X-direction first
        my_numx = nx/npcol       
        modcheck = MOD(nx,npcol)
        IF (my_row_rank .lt. modcheck) my_numx = my_numx+1
        my_imin = 1 + my_row_rank*my_numx
        IF (modcheck .le. my_row_rank) my_imin = my_imin+modcheck
        my_imax = my_imin+my_numx-1        

        ! Y-direction first
        my_numy = ny/nprow       
        modcheck = MOD(ny,nprow)
        IF (my_col_rank .lt. modcheck) my_numy = my_numy+1
        my_jmin = 1 + my_col_rank*my_numy
        IF (modcheck .le. my_col_rank) my_jmin = my_jmin+modcheck
        my_jmax = my_jmin+my_numy-1        


    END SUBROUTINE LOAD_BALANCE


    SUBROUTINE INIT_ARR(arr, amp, orderx, ordery, orderz)
        IMPLICIT NONE
        ! FORTRAN will assume arr is indexed beginning with 1 by default
        ! We can alter this for the x and y dimension in the declaration as follows
        REAL*8, INTENT(INOUT) :: arr(my_imin-1:,my_jmin-1:,:)
        REAL*8, INTENT(IN) :: amp
        INTEGER, INTENT(IN) :: orderx, ordery, orderz
        REAL*8 :: sinkx, sinky, sinkz
        REAL*8 :: kx, ky, kz
        REAL*8, PARAMETER :: pi = 3.1415926535897932384626433832795028841972d0
        INTEGER :: i,j,k, dims(3), ni,nj,nk
        dims = shape(arr)


        kx = orderx*(pi/(nx-1))
        ky = ordery*(pi/(ny-1))
        kz = orderz*(pi/(nz-1))

        DO k = 1, nz
            sinkz = 1.0d0 !sin(kz*(k-1)) We will expand to 3-D later
            DO j = my_jmin, my_jmax
                sinky = sin(ky*(j-1))
                DO i = my_imin,my_imax
                    sinkx = sin(kx*(i-1))
                    arr(i,j,k) = arr(i,j,k)+amp*sinkx*sinky*sinkz
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE INIT_ARR
END PROGRAM MAIN
