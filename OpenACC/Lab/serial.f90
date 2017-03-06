PROGRAM MAIN
    IMPLICIT NONE
    REAL*8, ALLOCATABLE :: vars(:,:,:), tmp(:,:,:)
    INTEGER ::  nx,ny,nz,nt

    INTEGER :: i,j,k,n
    INTEGER :: sf =2
    nx = 128*sf
    ny = 128*sf
    nz = 128*sf
    nt = 100

    ALLOCATE(vars(1:nx,1:ny,1:nz))
    ALLOCATE( tmp(1:nx,1:ny,1:nz))
    
    vars(:,:,:) = 0.0d0
    tmp(:,:,:)   = 0.0d0

    DO n = 1, nt,2
        Write(6,*)n
        CALL Laplacian(vars,tmp)
        Write(6,*)n+1
        CALL Laplacian(tmp,vars)
    ENDDO

    
CONTAINS

    SUBROUTINE Laplacian(arrin, arrout)
        IMPLICIT NONE
        REAL*8, INTENT(In) :: arrin(:,:,:)
        REAL*8, INTENT(Out) :: arrout(:,:,:)
        DO k = 2, nz-1
            DO j = 2, ny-1
                DO i = 2, nz-1
                    arrout(i,j,k) = arrin(i,j,k) + &
                            arrin(i-1,j,k) + arrin(i+1,j,k) + &
                            arrin(i,j-1,k) + arrin(i,j+1,k)  
                ENDDO
            ENDDO
        ENDDO
    END SUBROUTINE Laplacian

END PROGRAM MAIN
