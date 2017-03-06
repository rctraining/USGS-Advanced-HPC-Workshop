! Example 0:  Data is placed on the GPU and resides there for 
! the duration of the code's execution.
! We can think of this as an application designed for a single GPU.
! There is no MPI-related communication going on in the background.

PROGRAM MAIN
    IMPLICIT NONE
    REAL*8, ALLOCATABLE :: vars(:,:,:), tmp(:,:,:)
    INTEGER ::  nx,ny,nz,nt

    INTEGER :: i,j,k,n
    INTEGER :: sf = 2
    nx = 128*sf
    ny = 128*sf
    nz = 128*sf
    nt = 100

    ALLOCATE(vars(1:nx,1:ny,1:nz))
    ALLOCATE( tmp(1:nx,1:ny,1:nz))
    
    vars(:,:,:) = 0.0d0
     tmp(:,:,:) = 0.0d0
    !$acc data copy(vars, tmp)
    DO n = 1, nt,2
        WRITE(6,*)n
        CALL Laplacian(vars,tmp)
        WRITE(6,*)n+1
        CALL Laplacian(tmp,vars)
        

    ENDDO
    !$acc end data
    
CONTAINS

    SUBROUTINE Laplacian(arrin, arrout)
        IMPLICIT NONE
        REAL*8, INTENT(In) :: arrin(:,:,:)
        REAL*8, INTENT(Out) :: arrout(:,:,:)
        !$acc parallel loop present(arrin,arrout) collapse(3)
        DO k = 2, nz-1
            DO j = 2, ny-1
                DO i = 2, nz-1
                    arrout(i,j,k) = arrin(i,j,k) + &
                            arrin(i-1,j,k) + arrin(i+1,j,k) + &
                            arrin(i,j-1,k) + arrin(i,j+1,k)  
                ENDDO
            ENDDO
        ENDDO
        !$acc end parallel loop
    END SUBROUTINE Laplacian

END PROGRAM MAIN
