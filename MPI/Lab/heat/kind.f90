! To define a portable way of precision
! 
MODULE kind_module
! 
! Description: 
!   Define kinds for single and double precision in a portable way
! 
! Current Code Owner: Thomas Hauser 
! 
! History: 
!  
! Version   Date     Comment 
! -------   ----     ------- 
! 1.0       09-03-08 Original code. Thomas Hauser
! 
! Code Description: 
!   Language:           Fortran 90. 
!   Software Standards: "European Standards for Writing and  
!     Documenting Exchangeable Fortran 90 Code". 
! 
     
  IMPLICIT NONE 
! Global (i.e. public) Declarations: 
! Global Parameters: 
! Symbolic names for kind type of single and double-precision reals:
! (with at least 6 and 12 digits of accuracy)

  INTEGER, PARAMETER :: single = SELECTED_REAL_KIND(p=6,r=20)
  INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(p=12,r=150)
  INTEGER, PARAMETER :: prec = double

! Symbolic names for kind type of single and double-precision complex:

  INTEGER, PARAMETER :: cSingle = KIND((1.0_single,1.0_single))
  INTEGER, PARAMETER :: cDouble = KIND((1.0_double,1.0_double))
  INTEGER, PARAMETER :: cPrec = cSingle
         
END MODULE kind_module
         
