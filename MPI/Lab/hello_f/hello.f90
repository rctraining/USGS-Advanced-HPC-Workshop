PROGRAM hello


  INTEGER, PARAMETER:: dp=KIND(0.d0) 

  INTEGER ::  ierr, comm_size, comm_rank, length
  REAL(dp) :: RESULT
  CHARACTER (len=100) :: proc_name = "test"

  WRITE(*,*) "Hello World from process = ", comm_rank, " on processor ", proc_name

  RESULT = EXP(REAL(comm_rank, dp))

  WRITE(*,*) "Exp(", comm_rank, ") = ", RESULT 

  WRITE(*,*) "Number of mpi processes = ", comm_size

END PROGRAM hello
