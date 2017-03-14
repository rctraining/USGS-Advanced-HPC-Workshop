
PROGRAM hello


  IMPLICIT NONE

  INTEGER :: nthreads, tid

  nthreads = 1
  tid = 0
  WRITE(*,*) 'Hello World from thread = ', tid
  WRITE(*,*) 'Number of threads = ', nthreads

END PROGRAM hello
