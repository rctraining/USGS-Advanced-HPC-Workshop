PROGRAM MAIN
    IMPLICIT NONE
    REAL*8, ALLOCATABLE :: vars(:,:,:), tmp(:,:,:)
    INTEGER ::  nx,ny,nz,nt

    INTEGER :: i,j,k,n,nq
    INTEGER :: sf = 2
    nx = 128*sf
    ny = 128*sf
    nz = 128*sf
    nt = 100
    nq = 4

    ALLOCATE(vars(1:nx,1:ny,1:nz))
    ALLOCATE( tmp(1:nx,1:ny,1:nz))
    
    vars(:,:,:) = 0.0d0
     tmp(:,:,:) = 0.0d0
    !$ACC enter data copyin(vars, tmp)
    DO n = 1, nt
        

        CALL Laplacian(vars,tmp,nq)

        !//////////////////////////////////////////
        ! For parallel calculations, communication
        ! of some sort will need to be done outside
        ! of the central kernel(s)
        !
        ! We use the assignment operation below as
        ! a placeholder for this additional overhead
        !
        ! Here, we are forced to move data onto and off of the GPU.
        ! This version is very slow compared to the serial analog.
        vars = tmp
        WRITE(6,*)n

    ENDDO
    !$ACC exit data delete(vars,tmp)
    
CONTAINS

    SUBROUTINE Laplacian(arrin, arrout, numq)
        IMPLICIT NONE
        REAL*8, INTENT(In)  :: arrin(:,:,:)
        REAL*8, INTENT(Out) :: arrout(:,:,:)
        INTEGER, INTENT(In) :: numq
        INTEGER :: queue

        DO k = 2, nz-1
            queue = MOD(k,numq)+1
            !$ACC update device(arrin(1:4,1,k)) async(queue)
            !$acc parallel loop present(arrin,arrout) collapse(2) 
            DO j = 2, ny-1
                DO i = 2, nz-1
                    arrout(i,j,k) = arrin(i,j,k) + &
                            arrin(i-1,j,k) + arrin(i+1,j,k) + &
                            arrin(i,j-1,k) + arrin(i,j+1,k)  
                ENDDO
            ENDDO
            !$acc end parallel loop
            !$ACC update host(arrout(1:4,1,k)) async(queue)
        ENDDO
        !$ACC WAIT
    END SUBROUTINE Laplacian

END PROGRAM MAIN
