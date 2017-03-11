PROGRAM MAIN
    USE ISO_FORTRAN_ENV, ONLY : output_unit,real64
    USE MPI
    IMPLICIT NONE
    REAL*8, ALLOCATABLE :: var(:), tmp(:)
    INTEGER ::  nx,nt
    REAL( real64) :: ct_one, ct_two, ct_three
    REAL( real64) :: elapsed_time, loop_time, init_time
    INTEGER :: i,j,k,n, ierr

    !/////////////////////////////////////////
    ! These two variables are initialized in INIT_COMM below.
    INTEGER :: my_rank, ncpu  ! Rank of this process; total number of MPI processes
    INTEGER :: my_imin, my_imax ! the range of x-indices that this rank owns
    INTEGER :: my_numx  ! The number of-indices assinged to this processes
                         ! my_numx = my_imax-my_imin+1

    CALL cpu_time( time= ct_one)

    !////////////////////////////////////////////////////
    ! Set some default values for the problem parameters.
    nx = 4096   ! x-dimension of our 3-D array
    nt = 4096   ! number of time steps to integrate over


    !////////////////////////////////////////////////////
    ! Check to see if the user has overridden
    ! these defaults at the command line.
    ! Command line calling syntax:
    ! ./ex5.out -nx 256 -nt 100 
    CALL grab_args(nx,nt)

    CALL INIT_COMM()

    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,'(a)') ' ///////////////////////////'
        WRITE(output_unit,'(a)') '  1-D Diffusion Equation   '
        WRITE(output_unit,'(a)') '  --Problem Parameters-- '
        WRITE(output_unit,'(a,i0)')'   nx =', nx
        WRITE(output_unit,'(a,i0)')'   nt =', nt
        WRITE(output_unit,'(a,i0)')' ncpu =', ncpu
        WRITE(output_unit,'(a)') ' ///////////////////////////'
        WRITE(output_unit,'(a)') ' '
    ENDIF
    !/////////////////////////////////////////////////////
    ! Initialize var and the work array

    !NOTE: my_imin & my_numx are correct only for ncpu = 1
    !      It should be modified based on my_rank.
    my_numx = nx
    my_imin = 1
    my_imax = my_imin+my_numx-1

    DO i = 0, ncpu-1
        IF (i .eq. my_rank) THEN
            WRITE(6,*)'my_rank, my_imin, my_imax: ', my_rank, my_imin, my_imax
        ENDIF
        ! The use of a barrier here helps the line above to be printed out
        ! "in order"
        IF (ncpu .gt. 1) CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    ENDDO

    
    ALLOCATE(var(my_imin-1:my_imax+1))  ! We allow space for two
    ALLOCATE(tmp(my_imin:my_imax))  ! "ghost zones" in the var array
    
    var(:)   = 0.0d0
    tmp(:)   = 0.0d0
    CALL INIT_ARR(var, 1.0d0, 2) ! single sine wave in x-direction


    !/////////////////////////////////////////////////////
    ! Evolve the system
    CALL cpu_time( time= ct_two)
    DO n = 1, nt
        IF (my_rank .eq. 0) THEN
            IF (MOD(n,100) .eq. 0) Write(output_unit,'(a,i5)')' Timestep: ',n
        ENDIF

        !* EDIT THIS FUNCTION to communicate ghost zones
        !  to left and right neighbor ranks.
        IF (ncpu .gt. 1) CALL GHOST_ZONE_COMM()


        CALL Laplacian(var,tmp)

    ENDDO
    ! We may want to write out the values for debugging
    !Write(output_unit,*) var(:)   

    !//////////////////////////////////////////
    ! Write out some timing information.


    CALL cpu_time( time= ct_three)
    elapsed_time = ct_three-ct_one
    init_time = ct_two-ct_one
    loop_time = ct_three-ct_two
    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,'(a)')' Complete!'
        WRITE(output_unit,'( a, ES14.4, a)')'         Elapsed time: ', elapsed_time, ' seconds.'
        WRITE(output_unit, fmt= '( a, ES14.4, a)') '  Initialization time: ', init_time, ' seconds.'
        WRITE(output_unit, fmt= '( a, ES14.4, a)') '            Loop time: ', loop_time, ' seconds.'
    ENDIF

    ! Always call this at the end of the program so that all ranks exit together.
    ! Once your code is working, try removing this call.  The code will still work,
    ! but sometimes you may not see the timing information printed out because other
    ! ranks exit before rank 0 can write.
    CALL MPI_FINALIZE(ierr)
CONTAINS

    SUBROUTINE GHOST_ZONE_COMM()
        IMPLICIT NONE
        ! Communicate boundary values here...

    END SUBROUTINE GHOST_ZONE_COMM


    SUBROUTINE Laplacian(arrin, work)
        IMPLICIT NONE
        REAL*8, INTENT(InOut) :: arrin(my_imin-1:)
        REAL*8, INTENT(InOut) :: work(my_imin:)
        INTEGER :: ni, imin,imax

        ni = size(arrin)

        imin = max(my_imin,2)    ! Avoid the boundaries, which are held fixed at zero
        imax = min(my_imax,nx)

        DO i = imin, imax
            work(i) =  arrin(i-1) + arrin(i+1) 
        ENDDO


        DO i = imin, imax
            arrin(i) = work(i)*0.5d0
        ENDDO

    END SUBROUTINE Laplacian

    SUBROUTINE grab_args(numx, numiter)
            IMPLICIT NONE

            INTEGER, INTENT(OUT)   :: numx
            INTEGER, INTENT(OUT)   :: numiter


            INTEGER :: n                    ! Number of command-line arguments
            INTEGER :: i                    
            CHARACTER(len=1024) :: argname  ! Argument key
            CHARACTER(len=1024) :: val      ! Argument value



            n = command_argument_count()
            DO i=1,n,2
                    CALL get_command_argument(i, argname)
                    CALL get_command_argument(i+1, val)
                    SELECT CASE(argname)
                            CASE('-nx')
                                    read(val, '(I8)') numx
                            CASE('-nt')
                                    read(val, '(I8)') numiter
                            CASE DEFAULT
                                    WRITE(output_unit,'(a)') ' '
                                    WRITE(output_unit,'(a)') &
                                    ' Unrecognized option: '// trim(argname)
                    END SELECT
            ENDDO


    END SUBROUTINE grab_args

    SUBROUTINE INIT_COMM()
        IMPLICIT NONE
        INTEGER :: ierr
        !Initialize MPI  (ierr is an error flag with nonzero value if there's a problem)
        Call MPI_INIT( ierr )
        IF (ierr .ne. 0) STOP 'Error intializing MPI'

        !Find number of MPI processes executing this program (ncpu)
        Call MPI_Comm_size(MPI_COMM_WORLD, ncpu,ierr)
        IF (ierr .ne. 0) STOP 'Error finding ncpu'

        !Find this process's rank (my_rank; unique numeric identifier for each MPI process)
        Call MPI_Comm_rank(MPI_COMM_WORLD, my_rank,ierr)
        IF (ierr .ne. 0) STOP 'Error finding my_rank'
        WRITE(output_unit,'(a,i0,a,i0,a)')'Rank ',my_rank,' of ', ncpu, ' is active.'
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    END SUBROUTINE INIT_COMM


    SUBROUTINE INIT_ARR(arr, amp, orderx)
        IMPLICIT NONE
        REAL*8, INTENT(INOUT) :: arr(:)
        REAL*8, INTENT(IN) :: amp
        INTEGER, INTENT(IN) :: orderx
        REAL*8 :: sinkx
        REAL*8 :: kx
        REAL*8, PARAMETER :: pi = 3.1415926535897932384626433832795028841972d0
        INTEGER :: i,ni
        ni = size(arr)

        kx = orderx*(pi/(ni-1))

        DO i = my_imin, my_imax
            sinkx = sin(kx*(i-1))
            arr(i) = arr(i)+amp*sinkx
        ENDDO
        If (my_imin .eq. 1) arr(1) = 0
        If (my_imax .eq. ni) arr(ni) = 0

    END SUBROUTINE INIT_ARR
END PROGRAM MAIN
