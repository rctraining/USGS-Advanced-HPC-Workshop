/*!**********************************************************************
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
		      
int main(int argc, char *argv[])
{
  double pi25dt = 3.141592653589793238462643;

  double  pi, h, sum, x, f, a;
  int n, i, ierr;
  

  printf("Enter the number of intervals: \n");
  scanf("%d", &n);

    
  /* calculate the interval size */
  h = 1.0/ (double) n;

  sum  = 0.0;
  for (i = 1; i<= n; i++) {
    x = h * ((double) i - 0.5);
    sum += 4.0 / (1.0 + x*x);
  }
  pi = h * sum;

  printf("pi is approximately %.16f, Error is %.16f\n", pi, fabs(pi-pi25dt));
}


