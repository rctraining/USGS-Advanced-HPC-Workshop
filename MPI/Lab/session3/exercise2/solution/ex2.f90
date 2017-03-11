!///////////////////////////////////////
! Solution
! Calling syntax:  mpiexec -ncpu ./ex1.out -nprow X -npcol Y
! WHERE X*Y = ncpu
! Wave_amp, which sets the global maximum, is set in main_input

PROGRAM MAIN
    USE MPI
    USE ISO_FORTRAN_ENV, ONLY : output_unit
    IMPLICIT NONE
    INTEGER :: my_rank, ncpu ! global rank and number of MPI processes
    INTEGER :: my_row_rank, my_col_rank ! The rank within a row and column for each process respectively
    INTEGER :: nprow=2, npcol=2 ! Number of process rows and columns respectively
                             ! 0 <= my_row_rank <= npcol
                             ! 0 <= my_col_rank <= nprow
    
    INTEGER :: ierr, i, j
    INTEGER :: col_comm, row_comm ! subcommunicators for rows and columns
                                   ! row_comm connects all processes in the same row
                                   ! col_comm connects all processes in the same column
    CHARACTER*120 :: msg
    REAL*8 :: wave_amp= 1.0d0  ! Amplitude of the sine wave we initialize (read from main_input) 
    Namelist /Input_Namelist/ wave_amp

    ! We have some additional data related to our distributed array
    INTEGER, PARAMETER :: nx_local = 32, ny_local = 32
    INTEGER :: nx_global, ny_global  ! = nx_local*nprow;  = ny_local*npcol
    REAL*8  :: var(nx_local,ny_local)  
    REAL*8  :: local_max = 0.0d0, global_max = 0.0d0


    !/////////////////////// Initialization
    CALL GRAB_ARGS()
    CALL INIT_COMM()
    CALL READ_INPUT_DATA() ! rank 0 reads input and broadcasts

    nx_global = nx_local*nprow
    ny_global = ny_local*nprow
    var(:,:) = 0

    CALL INIT_ARR(var, wave_amp, 2, 2) ! single sine hump in each direction
    local_max = maxval(var)
    global_max = local_max

    !/////////////////////////////////////////////////////////////////
    !  Report on the process parameters
    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,'(a)') ' ///////////////////////////'
        WRITE(output_unit,'(a)') '  Row-reading program  '
        WRITE(output_unit,'(a)') '  --MPI Parameters-- '
        WRITE(output_unit,'(a,i0)')'   nprow = ', nprow
        WRITE(output_unit,'(a,i0)')'   npcol = ', npcol
        WRITE(output_unit,'(a,i0)')'   ncpu  = ', ncpu
        WRITE(output_unit,'(a)') ' ///////////////////////////'
        WRITE(output_unit,'(a)') ' '
        WRITE(output_unit,'(a)') ' Proccess Grid Layout: '
    ENDIF

    WRITE(msg,'(a,i0,a,i0,a,i0)') &
        ' row = ', my_col_rank, ' column = ', my_row_rank, ' rank = ', my_rank
    CALL ORDERED_PRINT(msg)

    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,*)' '
        WRITE(output_unit,*)'Max values prior to row_col_max...'
    ENDIF
    WRITE(msg,'(a,i0,a,i0,a,ES12.4,a,ES12.4)') &
        ' row = ', my_col_rank, ' column = ', my_row_rank, &
        ' local max = ', local_max, ' global max = ', global_max
    CALL ORDERED_PRINT(msg)

    CALL ROW_COL_MAX()

    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,*)' '
        WRITE(output_unit,*)'Max values following to row_col_max...'
    ENDIF

    WRITE(msg,'(a,i0,a,i0,a,ES12.4,a,ES12.4)') &
        ' row = ', my_col_rank, ' column = ', my_row_rank, &
        ' local max = ', local_max, ' global max = ', global_max
    CALL ORDERED_PRINT(msg)


    CALL MPI_FINALIZE(ierr)

CONTAINS

    SUBROUTINE ROW_COL_MAX()
        IMPLICIT NONE
        INTEGER :: ierr
        REAL*8 :: tmp
        ! Reduce local_max into tmp based on row-wise values
        CALL MPI_Allreduce(local_max, tmp, 1, MPI_DOUBLE_PRECISION, & 
            MPI_MAX, row_comm, ierr)
        ! Reduce tmp into global_max based on column-wise values
        CALL MPI_Allreduce(tmp, global_max, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, col_comm, ierr)  
    END SUBROUTINE ROW_COL_MAX

    SUBROUTINE ORDERED_PRINT(message)
        IMPLICIT NONE
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

    SUBROUTINE READ_INPUT_DATA()
        IMPLICIT NONE
        Character*120 :: input_file
        input_file = 'main_input'

		! Rank 0 reads from input
        IF (my_rank .eq. 0) THEN
		    OPEN(unit=20, file=input_file, status="old", position="rewind")
		    READ(unit=20, nml=input_namelist)
		    CLOSE(20)
        ENDIF

        ! Broadcast across column 0 and then across rows
        IF (my_row_rank .eq. 0) THEN
            CALL MPI_Bcast( wave_amp, 1, MPI_DOUBLE_PRECISION, 0,col_comm,ierr)
        ENDIF
        CALL MPI_Bcast( wave_amp, 1, MPI_DOUBLE_PRECISION, 0,row_comm,ierr) 
    END SUBROUTINE READ_INPUT_DATA

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
        INTEGER :: test
        my_row_rank = MOD(my_rank,npcol)
        my_col_rank = my_rank/npcol

        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, my_col_rank, my_rank, row_comm, ierr)
        !CALL mpi_comm_rank(row_comm, test, ierr)
        !WRITE(6,*)test, my_row_rank

        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, my_row_rank, my_rank, col_comm, ierr)
        !CALL mpi_comm_rank(col_comm, test, ierr)
        !WRITE(6,*)test, my_col_rank

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
                            CASE DEFAULT
                                    WRITE(output_unit,'(a)') ' '
                                    WRITE(output_unit,'(a)') &
                                    ' Unrecognized option: '// trim(argname)
                    END SELECT
            ENDDO


    END SUBROUTINE grab_args
    SUBROUTINE INIT_ARR(arr, amp, orderx, ordery)
        IMPLICIT NONE
        REAL*8, INTENT(INOUT) :: arr(:,:)
        REAL*8, INTENT(IN) :: amp
        INTEGER, INTENT(IN) :: orderx, ordery
        REAL*8 :: sinkx, sinky
        REAL*8 :: kx, ky
        REAL*8, PARAMETER :: pi = 3.1415926535897932384626433832795028841972d0
        INTEGER :: i,j,k, dims(2), ni,nj
        INTEGER :: iglobal, jglobal
        dims = shape(arr)
        nj = dims(2)
        ni = dims(1)

        kx = orderx*(pi/(nx_global-1))
        ky = ordery*(pi/(ny_global-1))
        
        DO j = 1, nj
            jglobal = my_col_rank*ny_local+j
            sinky = sin(ky*(jglobal-1))
            DO i = 1, ni
                iglobal = my_row_rank*nx_local+i
                sinkx = sin(kx*(iglobal-1))
                arr(i,j) = arr(i,j)+amp*sinkx*sinky
            ENDDO
        ENDDO

    END SUBROUTINE INIT_ARR


END PROGRAM MAIN
