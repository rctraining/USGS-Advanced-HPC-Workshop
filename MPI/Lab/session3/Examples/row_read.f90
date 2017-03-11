!///////////////////////////////////////
! Example program:   2-D process grid
! It is common practice to assume a virtual topology for a set of MPI ranks.
! In this example, we organize MPI ranks into a process grid of rows and columns.
!
! Processes in column 0 read from an input file and then broadcast the information
! to other processes within the row.  The input file is named row_read_input.
!
! I/O / Communication patterns such as this are one way to reduce the number of 
! processes accessing the filesystem...
!
! Calling syntax:  mpiexec -ncpu ./row_read.out -nprow X -npcol Y
! WHERE X*Y = ncpu

PROGRAM MAIN
    USE MPI
    USE ISO_FORTRAN_ENV, ONLY : output_unit
    IMPLICIT NONE
    INTEGER :: my_rank, ncpu ! global rank and number of MPI processes
    INTEGER :: my_row_rank, my_col_rank ! The rank within a row and column for each process respectively
    INTEGER :: nprow=1, npcol=1 ! Number of process rows and columns respectively
                             ! 0 <= my_row_rank <= npcol
                             ! 0 <= my_col_rank <= nprow
    
    INTEGER :: ierr, i, j
    INTEGER :: col_comm, row_comm ! subcommunicators for rows and columns
                                   ! row_comm connects all processes in the same row
                                   ! col_comm connects all processes in the same column
    REAL*8 :: input_data = 0.0d0
    Namelist /Input_Namelist/ input_data



    !/////////////////////// Initialization
    CALL GRAB_ARGS()
    CALL INIT_COMM()

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
    DO j = 0, nprow-1    ! Loop over each row
        DO i = 0, npcol-1 ! Loop over each column
            IF ((my_row_rank .eq. i) .and. (my_col_rank .eq. j) ) THEN
                WRITE(output_unit,'(a,i0,a,i0,a,i0)') &
                    ' row = ', my_col_rank, ' column = ', my_row_rank, ' rank = ', my_rank
            ENDIF
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ENDDO
    ENDDO

    !////////////////////////////////////////////
    !  READ INPUT and Broadcast
    IF (my_row_rank .eq. 0) THEN
        ! First process in each row reads the input file
        CALL READ_INPUT_DATA()
    ENDIF

    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,*)' '
        WRITE(output_unit,*)'Input_data prior to broadcast...'
    ENDIF

    DO j = 0, nprow-1    ! Loop over each row
        DO i = 0, npcol-1 ! Loop over each column
            IF ((my_row_rank .eq. i) .and. (my_col_rank .eq. j) ) THEN
                WRITE(output_unit,'(a,i0,a,i0,a,F8.4)') &
                    ' row = ', my_col_rank, ' column = ', my_row_rank, ' input_data = ', input_data
            ENDIF
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ENDDO
    ENDDO

    ! The value of input_data stored rank 0 of row_comm is broadcast
    ! to all other processes within the communicator.
    ! We could alternatively broadcast from rank 3, say, by changing the 0 to a 3 below...
    CALL MPI_Bcast( input_data, 1, MPI_DOUBLE_PRECISION, 0,row_comm,ierr) 

    IF (my_rank .eq. 0) THEN
        WRITE(output_unit,*)' '
        WRITE(output_unit,*)'Input_data following broadcast...'
    ENDIF

    DO j = 0, nprow-1    ! Loop over each row
        DO i = 0, npcol-1 ! Loop over each column
            IF ((my_row_rank .eq. i) .and. (my_col_rank .eq. j) ) THEN
                WRITE(output_unit,'(a,i0,a,i0,a,F8.4)') &
                    ' row = ', my_col_rank, ' column = ', my_row_rank, ' input_data = ', input_data
            ENDIF
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ENDDO
    ENDDO

    CALL MPI_FINALIZE(ierr)

CONTAINS

    SUBROUTINE READ_INPUT_DATA()
        IMPLICIT NONE
        Character*120 :: input_file
        input_file = 'row_read_input'

		! First read the main input file
		OPEN(unit=20, file=input_file, status="old", position="rewind")
		READ(unit=20, nml=input_namelist)
		CLOSE(20)

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

	!Function Init_SubGroup(ingrp,ncpus,ierr) result(grp)
        ! This routine is for multi-run jobs.  It splits ingrp into 
        ! an size(ncpus) subgroups, with the ith subgroup containing
        ! ncpus(i) ranks.
        !
        ! ingrp is a communicator with X ranks
        ! ncpus is an integer array arbitrary size <= X
        ! Sum(ncpus) must equal X (number of ranks in ingrp)

	!	Type(communicator) :: grp
!		Integer, Intent(out) :: ierr
    !    Integer, Intent(In) :: ncpus(1:)
    !    Type(communicator), Intent(In) :: ingrp
    !    Integer :: i,grank,gcolor, nclusters, mn_rank, mx_rank
    !    grank = ingrp%rank
    !    mn_rank = 0
    !    nclusters = size(ncpus)

    !    Do i = 1, nclusters
    !         mx_rank = ncpus(i)-1+mn_rank
    !         if ( (grank .ge. mn_rank) .and. (grank .le. mx_rank) ) Then
    !            gcolor = i-1
    !         Endif
    !         mn_rank = mx_rank+1
    !    Enddo
 !		Call mpi_comm_split(ingrp%comm, gcolor, ingrp%rank, grp%comm, ierr)!

    !	Call mpi_comm_size(grp%comm, grp%np, ierr)
	!	Call mpi_comm_rank(grp%comm, grp%rank, ierr)

	!End Function Init_SubGroup

END PROGRAM MAIN
