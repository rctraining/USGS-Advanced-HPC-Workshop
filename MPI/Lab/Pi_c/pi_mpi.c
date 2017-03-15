/*!**********************************************************************
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
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"
		      
int main(int argc, char *argv[])
{
  double pi25dt = 3.141592653589793238462643;

  double  mypi, pi, h, sum, x, f, a;
  int n, myid, numprocs, i, ierr;
  
  ierr = MPI_Init(&argc, &argv);
  
  ierr = MPI_Comm_size( MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myid);
  printf("Process %d of %d is running\n", myid, numprocs);

  if ( myid == 0 ) {
    printf("Enter the number of intervals: \n");
    scanf("%d", &n);
  }

  ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
  /* calculate the interval size */
  h = 1.0/ (double) n;

  sum  = 0.0;
  for (i = myid+1; i<= n; i += numprocs) {
    x = h * ((double) i - 0.5);
    sum += 4.0 / (1.0 + x*x);
  }
  mypi = h * sum;

  /* collect all the partial sums */
  ierr = MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* node 0 prints the answer. */
  if (myid == 0) {
    printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi-pi25dt));
  }

  ierr = MPI_Finalize();						     
}


