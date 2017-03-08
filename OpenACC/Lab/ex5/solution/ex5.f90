PROGRAM MAIN
    USE ISO_FORTRAN_ENV, ONLY : output_unit,real64
    USE OMP_LIB
    IMPLICIT NONE
    REAL*8, ALLOCATABLE :: var(:,:,:), tmp(:,:,:)
    INTEGER ::  nx,ny,nz,nt
    REAL( real64) :: ct_one, ct_two, ct_three
    REAL( real64) :: elapsed_time, loop_time, init_time
    INTEGER :: i,j,k,n, nq, nkgpu, nthread
    CALL cpu_time( time= ct_one)
    !////////////////////////////////////////////////////
    ! Set some default values for the problem parameters.
    nx = 256   ! x-dimension of our 3-D array
    ny = 256   ! y-dimension of our 3-D array
    nz = 256   ! z-dimension of our 3-D array
    nt = 100   ! number of time steps to integrate over
    nq = 4     ! The number of OpenACC queues to employ
    nkgpu = nz ! GPU handles z-levels from 2 through nkgpu-1
    !////////////////////////////////////////////////////
    ! Check to see if the user has specified anything 
    ! different at the command line.
    CALL grab_args(nx,ny,nz,nt,nq,nkgpu)

    ! IF the value of nk_GPU doesn't make sense, set it to nz.
    ! (so the GPU does all of the work)
    IF (nkGPU .lt. 2) nkGPU = nz
    IF (nkGPU .gt. nz) nkGPU = nz
    !$OMP PARALLEL
        nthread = omp_get_num_threads()
    !$OMP END PARALLEL

    WRITE(output_unit,'(a)') ' ///////////////////////////'
    WRITE(output_unit,'(a)') '  3-D Diffusion Equation   '
    WRITE(output_unit,'(a)') '  --Problem Parameters-- '
    WRITE(output_unit,'(a,i0)')'      nx =', nx
    WRITE(output_unit,'(a,i0)')'      ny =', ny
    WRITE(output_unit,'(a,i0)')'      nz =', nz
    WRITE(output_unit,'(a,i0)')'      nt =', nt
    WRITE(output_unit,'(a,i0)')'      nq =', nq
    WRITE(output_unit,'(a,i0)')'   nkgpu =', nkgpu 
    WRITE(output_unit,'(a,i0)')'   nthrd =', nthread
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


    !Write(6,*) var(nx/4-10:nx/4+10,ny/4,nz/4)

    Write(6,*) var(nx/4-2:nx/4+2,ny/4,5)
    Write(6,*) var(nx/4-2:nx/4+2,ny/4,nz-5)
    ! Create an initial copy of var and tmp on the GPU
    !$ACC enter data copyin(var, tmp)

    CALL cpu_time( time= ct_two)
    DO n = 1, nt
        IF (MOD(n,10) .eq. 0) Write(output_unit,'(a,i5)')' Timestep: ',n

        CALL Laplacian(var,tmp,nq,nkgpu)


        !/////////////////////////////////////
        ! This piece is carried out on the CPU 
        CALL GHOST_ZONE_COMM(var)
        
    ENDDO
    !Write(6,*) ' '
    !Write(6,*) var(nx/4-10:nx/4+10,ny/4,nz/4)
    !Write(6,*) ' '

    CALL cpu_time( time= ct_three)
    ! At the end, copy out the GPU portion of the var array and delete var & tmp on the GPU
    !$ACC update host(var(1:nx,1:ny,2:nkGPU-1))
    !$ACC exit data delete(var,tmp)

    Write(6,*) var(nx/4-2:nx/4+2,ny/4,5)
    Write(6,*) var(nx/4-2:nx/4+2,ny/4,nz-5)


    elapsed_time = ct_three-ct_one
    init_time = ct_two-ct_one
    loop_time = ct_three-ct_two
    WRITE(output_unit,'(a)')' Complete!'
    WRITE(output_unit,'( a, ES14.4, a)')'         Elapsed time: ', elapsed_time, ' seconds.'
    WRITE(output_unit, fmt= '( a, ES14.4, a)') '  Initialization time: ', init_time, ' seconds.'
    WRITE(output_unit, fmt= '( a, ES14.4, a)') '            Loop time: ', loop_time, ' seconds.'

    
CONTAINS

    SUBROUTINE Laplacian(arrin, work,numq,nk_GPU)
        IMPLICIT NONE
        REAL*8, INTENT(InOut) :: arrin(:,:,:)
        REAL*8, INTENT(InOut) :: work(:,:,:)
        Real*8 :: one_sixth
        INTEGER :: dims(3), ni, nj, nk, queue
        INTEGER, INTENT(In) :: numq, nk_GPU
        dims = shape(arrin)
        nk = dims(3)
        nj = dims(2)
        ni = dims(1)
        one_sixth = 1.0d0/6.0d0

        !////////////////////////////////////////////////////
        ! GPU portion of the derivative loop
        !$acc update device(arrin(1:nx,1:ny,1), arrin(1:nx,1:ny,nk_gpu))
        DO k = 2, nk_GPU-1
            queue = MOD(k,numq)+1

            !$acc  update device(arrin(1,1:ny,k),arrin(nx,1:ny,k), &
            !$acc    arrin(1:nx,1,k), arrin(1:nx,ny,k)) async(queue)
            !$acc parallel loop present(arrin,work) collapse(2) async(queue)
            DO j = 2, nj-1
                DO i = 2, ni-1
                    work(i,j,k) =  &
                            arrin(i-1,j,k) + arrin(i+1,j,k) + &
                            arrin(i,j-1,k) + arrin(i,j+1,k) + & 
                            arrin(i,j,k-1) + arrin(i,j,k+1)
                ENDDO
            ENDDO
            !$acc end parallel loop
        ENDDO

        !/////////////////////////////
        ! CPU portion of the derivative loop
        !$OMP PARALLEL SHARED(work,arrin), PRIVATE(i,j,k)
   
        !$OMP DO
        DO k = nk_GPU, nk-1
            DO j = 2, nj-1
                DO i = 2, ni-1
                    work(i,j,k) =  &
                            arrin(i-1,j,k) + arrin(i+1,j,k) + &
                            arrin(i,j-1,k) + arrin(i,j+1,k) + & 
                            arrin(i,j,k-1) + arrin(i,j,k+1)
                ENDDO
            ENDDO

        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL


        !///////////////////////////////////////////////////////
        ! Before updating arrin, we need to ensure that work has 
        ! been completely updated.  

        !$ACC WAIT

        !///////////////////////////////////////////////////////
        ! GPU portion of the update loop
        DO k = 2, nk_GPU-1
            queue = MOD(k,numq)+1
            !$acc parallel loop present(arrin,work) collapse(2) async(queue)
            DO j = 2, nj-1
                DO i = 2, ni-1
                    arrin(i,j,k) = work(i,j,k)*one_sixth
                ENDDO
            ENDDO
            !$acc end parallel loop
            !$acc  update host(arrin(2,1:ny,k),arrin(nx-1,1:ny,k), &
            !$acc    arrin(1:nx,2,k), arrin(1:nx,ny-1,k)) async(queue)
        ENDDO

        !////////////////////////////////////////////
        ! CPU Portion of the update loop
        !$OMP PARALLEL SHARED(work,arrin), PRIVATE(i,j,k)
        !$OMP DO
        DO k = nk_GPU, nk-1

            DO j = 2, nj-1
                DO i = 2, ni-1
                    arrin(i,j,k) = work(i,j,k)*one_sixth
                ENDDO
            ENDDO

        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL

        !//////////////////////////////////////////////////////////////////////////
        ! Wait until arrin has been updated on the GPU, then pull back over to CPU
        !$ACC WAIT

        !$acc update host(var(1:nx,1:ny,2), var(1:nx,1:ny,nk_GPU-1))
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


        DO k = 1, nk
            arr(1,:,k) = arr(ni-1,:,k)*0
            arr(:,1,k) = arr(:,nj-1,k)*0
            arr(ni,:,k) = arr(2,:,k)*0
            arr(:,nj,k) = arr(:,2,k)*0
        ENDDO
        arr(:,:,1) = arr(:,:,nk-1)*0
        arr(:,:,nk) = arr(:,:,2)*0

    END SUBROUTINE GHOST_ZONE_COMM


    SUBROUTINE grab_args(numx, numy, numz, numiter, numq, nk_gpu)
            IMPLICIT NONE

            INTEGER, INTENT(OUT)   :: numx
            INTEGER, INTENT(OUT)   :: numy
            INTEGER, INTENT(OUT)   :: numz
            INTEGER, INTENT(OUT)   :: numiter
            INTEGER, INTENT(INOUT) :: numq, nk_gpu


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
                            CASE('-nq')
                                    read(val, '(I8)') numq
                            CASE('-nkgpu')
                                    read(val, '(I8)') nk_gpu
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
            sinkz = sin(kz*k)
            DO j = 1, nj
                sinky = sin(ky*j)
                DO i = 1, ni
                    sinkx = sin(kx*i)
                    arr(i,j,k) = arr(i,j,k)+amp*sinkx*sinky*sinkz
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE INIT_ARR

END PROGRAM MAIN
