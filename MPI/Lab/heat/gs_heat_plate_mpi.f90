PROGRAM gs_heat_plate_mpi
  ! Description:
  !   This program solves the heat transfer in celsius of a plate in a square domain using the gauss-sidel method.
  ! Input files:
  !   The input to this file is inteded to come from stdin. However from the CLI a file
  !   is used with the < operator. In this file, the domain, heater_x, 
  !   and heater_y should be specified in that order on seperate lines for it work wo/ user input. 

  USE kind_module, ONLY : prec
  USE mympi
  USE mpi

  IMPLICIT NONE
  !
  ! Include statements
  ! Declarations must be of the form:
  ! [type] [VariableName]     !Description

  !Local parameters:

  !Local scalars:
  INTEGER :: i,j,iter=0 !- Loop variables
  REAL(kind = prec) :: elapsed_time
  INTEGER :: error
  LOGICAL :: quitProg = .FALSE.
  !
  ! The varaibles below are the initial temperatures for 
  ! the boudaries of the plate and the value of convergence.
  !
  TYPE heat_info_type
      INTEGER :: domain ! The area that contains the nodes of the virtual plate
      INTEGER :: heater_x, heater_y  ! The x and y location of the heater
      INTEGER :: heater_rank ! the processor rank that contains the heater
      REAL(kind = prec) :: init_temp !-initial temperature of the plate
      REAL(kind = prec) :: heater_temp !-initial temperature of the heater
      REAL(kind = prec) :: heat_top, heat_bottom, heat_left, heat_right!-initial temperatures of the boundaries
      REAL(kind = prec) :: convergence !-convergence criteria to end program
  END TYPE heat_info_type

  TYPE(heat_info_type) :: heat_info ! heat information for each processor

  TYPE grid_info_type
     INTEGER :: p    ! total number of processors
     INTEGER :: commcart  ! grid communicator
     INTEGER :: rank     ! processors rank in the grid
     INTEGER :: north    ! north neighbor
     INTEGER :: south   ! south neighbor
     INTEGER :: east     ! east neighbor
     INTEGER :: west    ! west neighbor
     INTEGER :: nx       ! number of local points in the x direction
     INTEGER :: ny       ! number of local points in the y direction
     ! Start of receive buffers
     REAL(kind = prec), ALLOCATABLE :: rnorth_buffer(:)
     REAL(kind = prec), ALLOCATABLE :: rsouth_buffer(:)
     REAL(kind = prec), ALLOCATABLE :: rwest_buffer(:)
     REAL(kind = prec), ALLOCATABLE :: reast_buffer(:)
     ! Start of send buffers
     REAL(kind = prec), ALLOCATABLE :: snorth_buffer(:)
     REAL(kind = prec), ALLOCATABLE :: ssouth_buffer(:)
     REAL(kind = prec), ALLOCATABLE :: swest_buffer(:)
     REAL(kind = prec), ALLOCATABLE :: seast_buffer(:)
  END TYPE grid_info_type

  ! Type instantiations
  TYPE(grid_info_type) :: grid_info ! grid information for each processor

  !- Local arrays:
  REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: t
  REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: t_old

  !
  ! MPI variables 
  !
  INTEGER :: ndummy = 1
  INTEGER :: globaliter = 0
  INTEGER :: mpierror, procs, myid !
  INTEGER :: nx !domain in x ;
  INTEGER :: ny !domain in y

  heat_info%domain = 500
  nx = heat_info%domain
  ny = heat_info%domain
!
!
! Initialize temperature variables
!
!
heat_info%heater_x = 251
heat_info%heater_y = 251
! check the bounds of the heaters location.
IF (heat_info%heater_x>heat_info%domain .OR. heat_info%heater_x<1 .OR. &
     &heat_info%heater_y < 1.OR. heat_info%heater_y>heat_info%domain) THEN
		WRITE(*,*) 'Error: invalid heater location - it is out of the domain.'
END IF

heat_info%init_temp = 21._prec !the initial temperature of the plate
heat_info%heater_temp = 29._prec !the initial temperature of the heater
!
! Start of boundary temperatures
!
heat_info%heat_top = 20._prec
heat_info%heat_bottom = 20._prec
heat_info%heat_left = 20._prec
heat_info%heat_right = 20._prec
!
!- End of boundary temperatures
!
heat_info%convergence = 0.0001_prec

  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(mpi_comm_world,procs,mpierror) ! store the number of processors
  CALL MPI_COMM_RANK(mpi_comm_world,myid,mpierror)

  ! print out the values used
  IF (myid == 0) THEN
	WRITE(*,*) "Values used: "
	WRITE(*,*) heat_info%domain, ' by', heat_info%domain, ' plate size'
	WRITE(*,*) 'Initial temperature of plate:', heat_info%init_temp
	WRITE(*,*) 'X location of heater:', heat_info%heater_x
	WRITE(*,*) 'Y location of heater:', heat_info%heater_y
	WRITE(*,*) 'Starting heater temperature:', heat_info%heater_temp
	WRITE(*,*) 'Top temperature:', heat_info%heat_top
	WRITE(*,*) 'Bottom temperature:', heat_info%heat_bottom
	WRITE(*,*) 'Right temperature:', heat_info%heat_right
	WRITE(*,*) 'Left temperature:', heat_info%heat_left
	WRITE(*,*) 'Convergence threshold used:',heat_info%convergence
        WRITE(*,*) 'Computing... Please wait... this could take a while'
  END IF

  elapsed_time = -mpi_wtime()

  CALL setup_grid(nx, ny, grid_info, heat_info)

  ! Allocate blocks on each processor
  ALLOCATE(t(1-ndummy:grid_info%nx+ndummy,1-ndummy:grid_info%ny+ndummy))
  ! Also make a second block to be able to check for convergence
  ALLOCATE(t_old(1-ndummy:grid_info%nx+ndummy,1-ndummy:grid_info%ny+ndummy))
  
  t(:,:) = REAL(heat_info%init_temp, prec) ! All the nodes have the same temperature
  
  ! setup heater 
  IF(grid_info%rank == heat_info%heater_rank) THEN
    t(heat_info%heater_x,heat_info%heater_y) = heat_info%heater_temp
  END IF
  
  run: DO
    t_old(:,:) = t(:,:)
    iter = iter+1
    !run one gauss-seidel iterations on the domain then go back to exchange bounds
    DO j = 1, grid_info%ny
       DO i = 1, grid_info%nx
          IF(grid_info%rank == heat_info%heater_rank) THEN
             IF(i == heat_info%heater_x .AND. j == heat_info%heater_y) CYCLE
          END IF
          t(i,j) = 0.25_prec*(t(i-1,j)+t(i+1,j)+t(i,j-1)+t(i,j+1))
       END DO
    END DO 
    CALL exchange_boundary(ndummy, grid_info,t,t_old, quitProg)
    ! check for convergence
    IF (quitProg .EQV. .TRUE.) THEN
	EXIT run
    END IF
  END DO run
  CALL mpi_barrier(grid_info%commcart,mpierror)

  elapsed_time = elapsed_time + mpi_wtime()

  CALL MPI_Allreduce(iter,globaliter,1,MPI_INTEGER, MPI_SUM, grid_info%commcart, mpierror)

  IF(grid_info%rank == 0) THEN
     WRITE(*,*)"Elapsed time",elapsed_time
     WRITE(*,*)"Number of iterations",globaliter
  ENDIF
  
  CALL mpi_finalize(mpierror)

CONTAINS

  SUBROUTINE setup_grid(gnx, gny, grid_info, heat_info)
   INTEGER, INTENT(in) :: gnx, gny ! size of domain in x and y
   TYPE(grid_info_type), INTENT(inout) :: grid_info
   TYPE(heat_info_type), INTENT(inout) :: heat_info
   INTEGER :: ndims = 2 ! number of dimensions in the grid
   LOGICAL :: reorder = .TRUE. ! do we want to re-order or not?
   INTEGER :: old_rank
   INTEGER :: psize(2) ! process grid dimensions
   INTEGER :: coords(2) !coordinates of the block in the grid
   LOGICAL :: periodic(2) = .FALSE. ! Is it perodic or not?
   INTEGER :: ierr
   CALL MPI_COMM_SIZE(mpi_comm_world,grid_info%p,ierr) ! store the number of processors
   CALL MPI_COMM_RANK(mpi_comm_world,old_rank,ierr) ! store the old ranks
   psize(:) = 0
   ! create a virtual grid
   CALL MPI_DIMS_CREATE(grid_info%p, ndims, psize, ierr)
   periodic(:) = .FALSE.
   ! create cartesian communicator        
   !call MPI_Cart_create(mpi_comm_world,ndims,psize,periodic,1,grid_info%commcart,ierr)
   CALL MPI_CART_CREATE(mpi_comm_world,ndims,psize,periodic,reorder,grid_info%commcart,ierr)
   !Get new rank
   CALL MPI_COMM_RANK(grid_info%commcart, grid_info%rank, ierr)
   CALL MPI_CART_COORDS(grid_info%commcart, grid_info%rank, ndims, coords, ierr) ! get the coords from the rank
   !write(*,*) "Rank",grid_info%rank,"coords",coords(1),coords(2)
   ! Build the blocks
   grid_info%nx = block_size(grid_info%rank, psize(1), gnx)
   grid_info%ny = block_size(grid_info%rank, psize(2), gny)
   ! Deterimine the processor with the heater
   coords(1) = (heat_info%heater_x)/(grid_info%nx)
   coords(2) = (heat_info%heater_y)/(grid_info%ny)
   CALL MPI_CART_RANK(grid_info%commcart, coords, heat_info%heater_rank, ierr)
   ! The actual heater value will have to be initialize after the block is created...
   ! Update heater_x and heater_y (the locations of the heater) with the offset needed to place the heater
   heat_info%heater_x = heat_info%heater_x - coords(1)*grid_info%nx
   heat_info%heater_y = heat_info%heater_y - coords(2)*grid_info%ny
   ! Setup north and south for communication
   CALL MPI_CART_SHIFT(grid_info%commcart, 0, 1, grid_info%north, grid_info%south,ierr)
   ! Setup east and west for communication
   CALL MPI_CART_SHIFT(grid_info%commcart, 1, 1, grid_info%west, grid_info%east, ierr)

   !allocate east_buffer
   ALLOCATE(grid_info%seast_buffer(grid_info%ny))
   ALLOCATE(grid_info%reast_buffer(grid_info%ny))
   !allocate west buffer
   ALLOCATE(grid_info%swest_buffer(grid_info%ny))
   ALLOCATE(grid_info%rwest_buffer(grid_info%ny))
   !allocate north buffer
   ALLOCATE(grid_info%snorth_buffer(grid_info%nx))
   ALLOCATE(grid_info%rnorth_buffer(grid_info%nx))
   !allocate south buffer
   ALLOCATE(grid_info%ssouth_buffer(grid_info%nx))
   ALLOCATE(grid_info%rsouth_buffer(grid_info%nx))

   !add heat information to dummy points which have no neighbors 
   IF(grid_info%north == MPI_PROC_NULL) THEN
      grid_info%rnorth_buffer=heat_info%heat_top
   ENDIF
   IF(grid_info%south == MPI_PROC_NULL) THEN
      grid_info%rsouth_buffer=heat_info%heat_bottom
   ENDIF
   IF(grid_info%east == MPI_PROC_NULL) THEN
      grid_info%reast_buffer=heat_info%heat_right
   ENDIF
   IF(grid_info%west == MPI_PROC_NULL) THEN
      grid_info%rwest_buffer=heat_info%heat_left
   ENDIF
   
END SUBROUTINE setup_grid

   SUBROUTINE exchange_boundary(ndummy, grid_info, t, t_old, quitProg)
     INTEGER, INTENT(in) :: ndummy
     TYPE(grid_info_type), INTENT(inout) :: grid_info
     REAL(kind = prec), INTENT(inout) :: t(1-ndummy:,1-ndummy:)
     REAL(kind = prec), INTENT(in) :: t_old(1-ndummy:,1-ndummy:)
     REAL(kind = prec) :: myDeltaT
     REAL(kind = prec) :: maxDeltaT 
     LOGICAL, INTENT(inout) :: quitProg
     INTEGER :: my_tag = 13
     INTEGER :: comrequest(8) ! 4 sends and 4 receives total of 8 communication operations 
     INTEGER :: ierr
     ! Setup a block to hold the previous value of t
     nx = grid_info%nx
     ny = grid_info%ny

     !initialize data that needs to be sent to other processors
     grid_info%seast_buffer = t(nx,1:ny)
     grid_info%swest_buffer = t(1,1:ny)
     grid_info%snorth_buffer = t(1:nx,1)
     grid_info%ssouth_buffer = t(1:nx,ny)
     ! generaly best to do recieves first
     CALL mpi_irecv(grid_info%reast_buffer, ny, MPI_DOUBLE_PRECISION, grid_info%east, my_tag,&
          & grid_info%commcart, comrequest(1), ierr)
     CALL mpi_irecv(grid_info%rwest_buffer, ny, MPI_DOUBLE_PRECISION, grid_info%west, my_tag, &
          &grid_info%commcart, comrequest(2), ierr)
     CALL mpi_irecv(grid_info%rnorth_buffer, nx, MPI_DOUBLE_PRECISION, grid_info%north, my_tag,&
          & grid_info%commcart, comrequest(3), ierr)
     CALL mpi_irecv(grid_info%rsouth_buffer, nx, MPI_DOUBLE_PRECISION, grid_info%south, my_tag,&
          & grid_info%commcart, comrequest(4), ierr)

     ! Now do all the sending
     CALL mpi_isend(grid_info%seast_buffer, ny, MPI_DOUBLE_PRECISION, grid_info%east, my_tag,&
          & grid_info%commcart, comrequest(5), ierr)
     CALL mpi_isend(grid_info%swest_buffer, ny, MPI_DOUBLE_PRECISION, grid_info%west, my_tag,&
          & grid_info%commcart, comrequest(6), ierr)
     CALL mpi_isend(grid_info%snorth_buffer, nx, MPI_DOUBLE_PRECISION, grid_info%north, my_tag,&
          & grid_info%commcart, comrequest(7), ierr)
     CALL mpi_isend(grid_info%ssouth_buffer, nx, MPI_DOUBLE_PRECISION, grid_info%south, my_tag,&
          & grid_info%commcart, comrequest(8), ierr)

     ! do something here to hide network speed :) Maybe check for convergance?
     ! each block will compute a finite differnce
     myDeltaT = MAXVAL(ABS(t(1:grid_info%nx,1:grid_info%ny) - t_old(1:grid_info%nx,1:grid_info%ny)))
     ! find the max difference value
     CALL MPI_Allreduce(myDeltaT,maxDeltaT,1,MPI_DOUBLE_PRECISION, MPI_MAX, grid_info%commcart, ierr)
     ! When the max difference is less than the convergance the program needs to quit
     WRITE (*,*) "Change in temperature: ", maxDeltaT
     IF (maxDeltaT < heat_info%convergence) THEN
        quitProg = .TRUE.
     END IF

     CALL mpi_waitall(8, comrequest, MPI_STATUSES_IGNORE, ierr) !this will become a dead lock if one of the processors is dead

     ! setup dummmy points with buffer information
     t(nx+1,1:ny)=grid_info%reast_buffer
     t(0,1:ny)=grid_info%rwest_buffer
     t(1:nx,ny+1)=grid_info%rsouth_buffer
     t(1:nx,0)=grid_info%rnorth_buffer
   END SUBROUTINE exchange_boundary
END PROGRAM gs_heat_plate_mpi
