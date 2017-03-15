!**********************************************************************
!   pi3f90.f - compute pi by integrating f(x) = 4/(1 + x**2)     
!
!  (C) 2001 by Argonne National Laboratory.
!  modified by Thomas Hauser
!     
!   Each node: 
!    1) receives the number of rectangles used in the approximation.
!    2) calculates the areas of it's rectangles.
!    3) Synchronizes for a global summation.
!   Node 0 prints the result.
!
!  Variables:
!
!    pi  the calculated result
!    n   number of points of integration.  
!    x           midpoint of each rectangle's interval
!    f           function to integrate
!    sum,pi      area of rectangles
!    tmp         temporary scratch space for global summation
!    i           do loop index
!****************************************************************************
PROGRAM pi_mpi

  USE mpi

  INTEGER, parameter :: dp = kind(1.0d0)
  REAL(dp)::  PI25DT = 3.141592653589793238462643d0

  REAL(dp) ::  mypi, pi, h, sum, x, f, a
  INTEGER :: n, myid, numprocs, i, rc
  !                                 function to integrate
  f(a) = 4.d0 / (1.d0 + a*a)

  CALL MPI_INIT( ierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
  PRINT *, 'Process ', myid, ' of ', numprocs, ' is alive'

  sizetype   = 1
  sumtype    = 2

  IF ( myid .EQ. 0 ) THEN
     WRITE(6,98)
98   FORMAT('Enter the number of intervals:')
     READ(5,99) n
99   FORMAT(i10)
  ENDIF

  CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  !                                 check for quit signal

  !                                 calculate the interval size
  h = 1.0d0/n

  sum  = 0.0d0
  DO i = myid+1, n, numprocs
     x = h * (DBLE(i) - 0.5d0)
     sum = sum + f(x)
  ENDDO
  mypi = h * sum

  !                                 collect all the partial sums
  CALL MPI_REDUCE(mypi,pi,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
       MPI_COMM_WORLD,ierr)

  !                                 node 0 prints the answer.
  IF (myid .EQ. 0) THEN
     WRITE(6, 97) pi, ABS(pi - PI25DT)
97   FORMAT('  pi is approximately: ', F18.16, &
          '  Error is: ', F18.16)
  ENDIF


  CALL MPI_FINALIZE(rc)
END PROGRAM pi_mpi




