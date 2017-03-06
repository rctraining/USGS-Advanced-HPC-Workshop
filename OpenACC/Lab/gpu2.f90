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
    !$ACC enter data copyin(vars, tmp)
    DO n = 1, nt
        
        !$ACC update device(vars)
        CALL Laplacian(vars,tmp)
        !$ACC update host(tmp)
        !//////////////////////////////////////////
        ! For parallel calculations, communication
        ! of some sort will need to be done outside
        ! of the central kernel(s)
        !
        ! We use the assignment operation below as
        ! a placeholder for this additional overhead
        !
        ! Here, we are forced to move data onto and off of the GPU.

        ! In this example, we only copy vars TO the GPU
        ! and we only copy tmp FROM the GPU.

        ! This version is still slow compared to the serial analog.
        vars = tmp
        WRITE(6,*)n

    ENDDO
    !$ACC exit data delete(vars,tmp)
    
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
