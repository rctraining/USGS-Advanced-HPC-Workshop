!******************************************************************************
! OpenMP Example - Hello World - Fortran Version
! FILE: hello.f90
! DESCRIPTION:
!   In this simple example, the master thread forks a parallel region.
!   All threads in the team obtain their unique thread number and print it.
!   The master thread only prints the total number of threads.  Two OpenMP
!   library routines are used to obtain the number of threads and each
!   thread's number.
! SOURCE: Blaise Barney  5/99
! LAST REVISED: Thomas Hauser 09/09/04
!******************************************************************************

PROGRAM hello

  use OMP_LIB , only : omp_get_thread_num, omp_get_num_threads

  IMPLICIT NONE

  INTEGER :: nthreads, tid

  !$OMP PARALLEL private(tid, nthreads)
  tid = omp_get_thread_num()
  WRITE(*,*) 'Hello World from thread = ', tid
  !$OMP BARRIER
  IF (tid == 0) THEN
     nthreads = omp_get_num_threads()
     WRITE(*,*) 'Number of threads = ', nthreads
  END IF
  !$OMP END PARALLEL

END PROGRAM hello
