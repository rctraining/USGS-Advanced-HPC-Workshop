PROGRAM MAIN
    USE ISO_FORTRAN_ENV, ONLY : output_unit,real64
    IMPLICIT NONE
    REAL*8, ALLOCATABLE :: var(:,:,:), tmp(:,:,:)
    INTEGER ::  nx,ny,nz,nt
    REAL( real64) :: ct_one, ct_two, ct_three
    REAL( real64) :: elapsed_time, loop_time, init_time
    INTEGER :: i,j,k,n
    CALL cpu_time( time= ct_one)
    !////////////////////////////////////////////////////
    ! Set some default values for the problem parameters.
    nx = 256   ! x-dimension of our 3-D array
    ny = 256   ! y-dimension of our 3-D array
    nz = 256   ! z-dimension of our 3-D array
    nt = 100   ! number of time steps to integrate over

    !////////////////////////////////////////////////////
    ! Check to see if the user has specified anything 
    ! different at the command line.

    ! Command line calling syntax:
    ! ./ex5.gpu -nx 256 -ny 256 -nz 256 -nt 100
    CALL grab_args(nx,ny,nz,nt)

    WRITE(output_unit,'(a)') ' ///////////////////////////'
    WRITE(output_unit,'(a)') '  3-D Diffusion Equation   '
    WRITE(output_unit,'(a)') '  --Problem Parameters-- '
    WRITE(output_unit,'(a,i0)')'   nx =', nx
    WRITE(output_unit,'(a,i0)')'   ny =', ny
    WRITE(output_unit,'(a,i0)')'   nz =', nz
    WRITE(output_unit,'(a,i0)')'   nt =', nt
    WRITE(output_unit,'(a)') ' ///////////////////////////'
    WRITE(output_unit,'(a)') ' '
    !/////////////////////////////////////////////////////
    ! Initialize var and the work array
    ALLOCATE(var(1:nx,1:ny,1:nz))
    ALLOCATE(tmp(1:nx,1:ny,1:nz))
    
    var(:,:,:) = 0.0d0
    tmp(:,:,:)   = 0.0d0
    CALL INIT_ARR(var, 1.0d0, 2, 2, 2) ! single sine wave in each dimension
    !/////////////////////////////////////////////////////
    ! Evolve the system:
    ! When running a parallel code, we inevitably need to pull data
    ! off of the GPU and carry out communication or non-GPU operations.
    ! The GHOST_ZONE_COMM function below serves as a placeholder
    ! for this communication and is assumed to be carried out on
    ! the CPU.  Modify the data movement portions of the code to
    ! reflect the fact that var should be updated on the CPU prior
    ! to ghost-zone communication.




    !$acc data copy(var, tmp)
    CALL cpu_time( time= ct_two)
    DO n = 1, nt

        IF (MOD(n,10) .eq. 0) Write(output_unit,'(a,i5)')' Timestep: ',n
        CALL Laplacian(var,tmp)

        !/////////////////////////////////////
        ! This piece is carried out on the CPU 
        CALL GHOST_ZONE_COMM(var)
        
    ENDDO
    CALL cpu_time( time= ct_three)
    !$acc end data

    Write(output_unit,*) var(8,8,8), var(nx-8,ny-8,nz-8)

    elapsed_time = ct_three-ct_one
    init_time = ct_two-ct_one
    loop_time = ct_three-ct_two
    WRITE(output_unit,'(a)')' Complete!'
    WRITE(output_unit,'( a, ES14.4, a)')'         Elapsed time: ', elapsed_time, ' seconds.'
    WRITE(output_unit, fmt= '( a, ES14.4, a)') '  Initialization time: ', init_time, ' seconds.'
    WRITE(output_unit, fmt= '( a, ES14.4, a)') '            Loop time: ', loop_time, ' seconds.'

    
CONTAINS

    SUBROUTINE Laplacian(arrin, work)
        IMPLICIT NONE
        REAL*8, INTENT(InOut) :: arrin(:,:,:)
        REAL*8, INTENT(InOut) :: work(:,:,:)
        Real*8 :: one_sixth
        INTEGER :: dims(3), ni, nj, nk
        dims = shape(arrin)
        nk = dims(3)
        nj = dims(2)
        ni = dims(1)
        one_sixth = 1.0d0/6.0d0
        !$acc parallel loop present(arrin,work) collapse(3)
        DO k = 2, nk-1
            DO j = 2, nj-1
                DO i = 2, ni-1
                    work(i,j,k) =  &
                            arrin(i-1,j,k) + arrin(i+1,j,k) + &
                            arrin(i,j-1,k) + arrin(i,j+1,k) + & 
                            arrin(i,j,k-1) + arrin(i,j,k+1)
                ENDDO
            ENDDO
        ENDDO
        !$acc end parallel loop

        !$acc parallel loop present(arrin,work) collapse(3)
        DO k = 2, nk-1
            DO j = 2, nj-1
                DO i = 2, ni-1
                    arrin(i,j,k) = work(i,j,k)*one_sixth
                ENDDO
            ENDDO
        ENDDO
        !$acc end parallel loop

    END SUBROUTINE Laplacian

    SUBROUTINE GHOST_ZONE_COMM(arr)
        IMPLICIT NONE
        REAL*8, INTENT(INOUT) :: arr(:,:,:)
        INTEGER :: dims(3)
        INTEGER :: i, j, k
        INTEGER :: ni, nj, nk
        dims = shape(arr)
        ni = dims(1)
        nj = dims(2)
        nk = dims(3)
        ! This is where we would normally communicate boundary
        ! information.
        ! Rather than invoking MPI for this example, 
        ! we simply copy the rightmost boundary to the leftmost boundary...


        DO k = 2, nk-1
            arr(1,:,k) = arr(ni-1,:,k)*0
            arr(:,1,k) = arr(:,nj-1,k)*0
            arr(ni,:,k) = arr(2,:,k)*0
            arr(:,nj,k) = arr(:,2,k)*0
        ENDDO
        arr(:,:,1) = arr(:,:,nk-1)*0
        arr(:,:,nk) = arr(:,:,2)*0

    END SUBROUTINE GHOST_ZONE_COMM




    SUBROUTINE grab_args(numx, numy, numz, numiter)
            IMPLICIT NONE

            INTEGER, INTENT(OUT)   :: numx
            INTEGER, INTENT(OUT)   :: numy
            INTEGER, INTENT(OUT)   :: numz
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
                            CASE('-ny')
                                    read(val, '(I8)') numy
                            CASE('-nz')
                                    read(val, '(I8)') numz
                            CASE('-nt')
                                    read(val, '(I8)') numiter
                            CASE DEFAULT
                                    WRITE(output_unit,'(a)') ' '
                                    WRITE(output_unit,'(a)') &
                                    ' Unrecognized option: '// trim(argname)
                    END SELECT
            ENDDO

    END SUBROUTINE grab_args

    SUBROUTINE INIT_ARR(arr, amp, orderx, ordery, orderz)
        IMPLICIT NONE
        REAL*8, INTENT(INOUT) :: arr(:,:,:)
        REAL*8, INTENT(IN) :: amp
        INTEGER, INTENT(IN) :: orderx, ordery, orderz
        REAL*8 :: sinkx, sinky, sinkz
        REAL*8 :: kx, ky, kz
        REAL*8, PARAMETER :: pi = 3.1415926535897932384626433832795028841972d0
        INTEGER :: i,j,k, dims(3), ni,nj,nk
        dims = shape(arr)
        nk = dims(3)
        nj = dims(2)
        ni = dims(1)

        kx = orderx*(pi/(ni-1))
        ky = ordery*(pi/(nj-1))
        kz = orderz*(pi/(nk-1))

        DO k = 1, nk
            sinkz = sin(kz*(k-1))
            DO j = 1, nj
                sinky = sin(ky*(j-1))
                DO i = 1, ni
                    sinkx = sin(kx*(i-1))
                    arr(i,j,k) = arr(i,j,k)+amp*sinkx*sinky*sinkz
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE INIT_ARR
END PROGRAM MAIN
