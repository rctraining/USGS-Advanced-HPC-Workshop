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

    !NOTE:  This is a worst-case example (motivated by my experience with 
    ! spectral/incompressible fluid solvers).  
    ! gpu3 copies data on/off, but only ghost zones.
    ! That's probably the better way to go.
    ! Maybe though, have them copy everything here at first to see why
    ! it's a bad idea...
    DO n = 1, nt
        !$acc data copy(vars, tmp)
        CALL Laplacian(vars,tmp)
        !$acc end data
        !//////////////////////////////////////////
        ! For parallel calculations, communication
        ! of some sort will need to be done outside
        ! of the central kernel(s)
        !
        ! We use the assignment operation below as
        ! a placeholder for this additional overhead
        !
        ! Here, we are forced to move data onto and off of the GPU
        ! during each iteration.   As a result, this version
        ! is very slow compared to the serial analog.
        vars = tmp
        WRITE(6,*)n

    ENDDO

    
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
