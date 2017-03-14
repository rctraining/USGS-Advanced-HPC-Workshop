!**********************************************************************
!   pi - compute pi by integrating f(x) = 4/(1 + x**2)     
!
!  (C) 2001 by Argonne National Laboratory. MPICH example
!  Modified by Thomas Hauser
!     
!   Program: 
!    1) reads the number of rectangles used in the approximation.
!    2) calculates the areas of it's rectangles.
!   
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
PROGRAM approx_pi

  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  REAL(dp)::  PI25DT = 3.141592653589793238462643d0
  
  REAL(dp) ::  pi, h, sum, x, f, a
  INTEGER :: n, i
  !                                 function to integrate
  f(a) = 4.0_dp / (1.0_dp + a*a)

  WRITE(6, 98)
98 FORMAT('Enter the number of intervals:')
  READ(5,99) n
99 FORMAT(i10)


  ! calculate the interval size
  h = 1.0_dp/ DBLE(n)

  sum  = 0.0_DP
  !$OMP PARALLEL DO PRIVATE(x), SHARED(h,n, sum)
  DO i = 1, n
     x = h * (DBLE(i) - 0.5_dp)
     !$OMP CRITICAL
     sum = sum + f(x)
     !$OMP END CRITICAL
  ENDDO
  !$OMP END PARALLEL DO
  pi = h * sum


  WRITE(6, 97) n, pi, ABS(pi - PI25DT)
97 FORMAT('  pi(', I15 ,' ) is approximately: ', F18.16, &
       '  Error is: ', F18.16)

END PROGRAM approx_pi




