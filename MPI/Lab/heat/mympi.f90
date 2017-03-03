MODULE mympi
  ! 
  ! Description: 
  !   Module for mpi to help associating data to multiple processors 
  ! 
  ! Current Code Owner: [Name of person responsible for this code] 
  ! 
  ! History: 
  !  
  ! Version   Date     Comment 
  ! -------   ----     ------- 
  ! 1.0       11-11-08 Original code. Thomas Hauser
  ! 
  ! Code Description: 
  !   Language:           Fortran 90. 
  !   Software Standards: "European Standards for Writing and  
  !     Documenting Exchangeable Fortran 90 Code". 
  ! 

  IMPLICIT NONE 

CONTAINS 

  INTEGER FUNCTION block_owner(j,p,n)
    INTEGER, INTENT(in) :: j,p,n
    block_owner = ((p*(j+1)-1)/n)
    RETURN
  END FUNCTION block_owner

  INTEGER FUNCTION block_low(id,p,n)
    INTEGER, INTENT(in) :: id,p,n
    block_low =  (id*n)/p       
    RETURN
  END FUNCTION block_low

  INTEGER FUNCTION block_high(id,p,n)
    INTEGER, INTENT(in) :: id,p,n
    block_high = ((id+1)*n)/p-1
    RETURN
  END FUNCTION block_high

  INTEGER FUNCTION block_size(id,p,n)
    INTEGER, INTENT(in) :: id,p,n
    block_size = (block_low(id+1,p,n)-block_low(id,p,n))
    RETURN
  END FUNCTION block_size

END MODULE mympi

