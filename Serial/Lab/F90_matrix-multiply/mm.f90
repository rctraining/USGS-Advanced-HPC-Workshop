!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! matrix_multiply.f90
!! Base matrix multiply code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM matrix_multiply
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: A(:,:),B(:,:),C(:,:)
  DOUBLE PRECISION :: FLOPS, MFLOPS,WT
  INTEGER :: n1, n2, n3
  character(len=100) :: myfmt

  write(*,*) "Matrix multiply: Fortran 90 version"
  write(*,*) "Compute matrix product C(n1, n3) = A(n1, n2) * B(n2, n3)"
  write(*,*) "Enter n1"
  read(*,*) n1
  write(*,*) "Enter n2"
  read(*,*) n2
  write(*,*) "Enter n3"
  read(*,*) n3

  myfmt = '(A, I7, A, I7, A)'
  write(*, myfmt) "Matrix A(", n1, ", ", n2, ")"
  write(*, myfmt) "Matrix B(", n2, ", ", n3, ")"
  write(*, myfmt) "Matrix C(", n1, ", ", n3, ")"

  flops = DBLE(n1) * DBLE(n2) * DBLE(n3)
  myfmt = '(A, ES15.2)'
  write(*, myfmt) "Number of floating point operations: ", flops


  ! C(n1, n3) = A(n1, n2) * B(n2, n3)
  ALLOCATE(A(n1, n2), B(n2, n3), C(n1, n3))

  call init_matrix(n1, n2, A)
  call init_matrix(n2, n3, B)

  ! Test for correct matrix dimensions
  write(*,*) "Method     Walltime [s]     MFLOPS        "
  write(*,*) "------   --------------  ----------------"
  myfmt = '(A6, F15.2, F15.2)'
  CALL ijk_mat_mult(n1, n2, n3, A, B, C, WT)
  MFLOPS = 2.d0*flops/(WT*1.d6)
  write(*, myfmt) 'ijk', wt, mflops

  CALL ikj_mat_mult(n1, n2, n3, A, B, C, WT)
  MFLOPS = 2.d0*flops/(WT*1.d6)
  write(*, myfmt) 'ikj', wt, mflops

  CALL jik_mat_mult(n1, n2, n3, A, B, C, WT)
  MFLOPS = 2.d0*flops/(WT*1.d6)
  write(*, myfmt) 'jik', wt, mflops

  CALL kij_mat_mult(n1, n2, n3, A, B, C, WT)
  MFLOPS = 2.d0*flops/(WT*1.d6)
  write(*, myfmt) 'kij', wt, mflops

  CALL kji_mat_mult(n1, n2, n3, A, B, C, WT)
  MFLOPS = 2.d0*flops/(WT*1.d6)
  write(*, myfmt) 'kji', wt, mflops

  CALL jki_mat_mult(n1, n2, n3, A, B, C, WT)
  MFLOPS = 2.d0*flops/(WT*1.d6)
  write(*, myfmt) 'jki', wt, mflops

  CALL f90_mat_mult(n1, n2, n3, A, B, C, WT)
  MFLOPS = 2.d0*flops/(WT*1.d6)
  write(*, myfmt) 'f90', wt, mflops

  DEALLOCATE(A,B,C)

contains

  subroutine init_matrix(n1, n2, m)
    integer, intent(in) :: n1, n2
    double precision, intent(out) :: m(n1, n2)

    integer, parameter :: s = 12345
    integer :: seed_size
    integer, allocatable :: seed(:)

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = s
    call random_seed(put=seed)
    call random_number(m)

    deallocate(seed)
  end subroutine init_matrix

  SUBROUTINE ijk_mat_mult(n1, n2, n3, A, B, C, WT)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, n2, n3
    DOUBLE PRECISION, DIMENSION(n1, n2), INTENT(in) :: A
    DOUBLE PRECISION, DIMENSION(n2, n3), INTENT(in) :: B
    DOUBLE PRECISION, DIMENSION(n1, n3), INTENT(out) :: C
    DOUBLE PRECISION, INTENT(out) :: WT
    double precision :: start_time, end_time
    INTEGER :: i, j, k

    C(:,:) = 0.0

    CALL cpu_time(start_time)

    DO i = 1, n1
       do j = 1, n3
          DO k = 1, n2
             C(i,j) = C(i,j) + A(i,k) * B(k,j)
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_time(end_time)

    WT = end_time - start_time

  END SUBROUTINE ijk_mat_mult


  SUBROUTINE ikj_mat_mult(n1, n2, n3, A, B, C, WT)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, n2, n3
    DOUBLE PRECISION, DIMENSION(n1, n2), INTENT(in) :: A
    DOUBLE PRECISION, DIMENSION(n2, n3), INTENT(in) :: B
    DOUBLE PRECISION, DIMENSION(n1, n3), INTENT(out) :: C
    DOUBLE PRECISION, INTENT(out) :: WT
    double precision :: start_time, end_time
    INTEGER :: i, j, k

    C(:,:) = 0.0

    CALL cpu_time(start_time)

    DO i = 1, n1
       DO k = 1, n2
          do j = 1, n3
             C(i,j) = C(i,j) + A(i,k) * B(k,j)
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_time(end_time)

    WT = end_time - start_time

  END SUBROUTINE ikj_mat_mult


  SUBROUTINE jik_mat_mult(n1, n2, n3, A, B, C, WT)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, n2, n3
    DOUBLE PRECISION, DIMENSION(n1, n2), INTENT(in) :: A
    DOUBLE PRECISION, DIMENSION(n2, n3), INTENT(in) :: B
    DOUBLE PRECISION, DIMENSION(n1, n3), INTENT(out) :: C
    DOUBLE PRECISION, INTENT(out) :: WT
    double precision :: start_time, end_time
    INTEGER :: i, j, k

    C(:,:) = 0.0

    CALL cpu_time(start_time)

    do j = 1, n3
       DO i = 1, n1
          DO k = 1, n2
             C(i,j) = C(i,j) + A(i,k) * B(k,j)
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_time(end_time)

    WT = end_time - start_time

  END SUBROUTINE jik_mat_mult

  SUBROUTINE kij_mat_mult(n1, n2, n3, A, B, C, WT)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, n2, n3
    DOUBLE PRECISION, DIMENSION(n1, n2), INTENT(in) :: A
    DOUBLE PRECISION, DIMENSION(n2, n3), INTENT(in) :: B
    DOUBLE PRECISION, DIMENSION(n1, n3), INTENT(out) :: C
    DOUBLE PRECISION, INTENT(out) :: WT
    double precision :: start_time, end_time
    INTEGER :: i, j, k

    C(:,:) = 0.0

    CALL cpu_time(start_time)

    DO k = 1, n2
       DO i = 1, n1
          do j = 1, n3
             C(i,j) = C(i,j) + A(i,k) * B(k,j)
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_time(end_time)

    WT = end_time - start_time

  END SUBROUTINE kij_mat_mult


  SUBROUTINE kji_mat_mult(n1, n2, n3, A, B, C, WT)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, n2, n3
    DOUBLE PRECISION, DIMENSION(n1, n2), INTENT(in) :: A
    DOUBLE PRECISION, DIMENSION(n2, n3), INTENT(in) :: B
    DOUBLE PRECISION, DIMENSION(n1, n3), INTENT(out) :: C
    DOUBLE PRECISION, INTENT(out) :: WT
    double precision :: start_time, end_time
    INTEGER :: i, j, k

    C(:,:) = 0.0

    CALL cpu_time(start_time)

    DO k = 1, n2
       DO j = 1, n3
          DO i = 1, n1
             C(i,j) = C(i,j) + A(i,k) * B(k,j)
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_time(end_time)

    WT = end_time - start_time

  END SUBROUTINE kji_mat_mult

  SUBROUTINE jki_mat_mult(n1, n2, n3, A, B, C, WT)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, n2, n3
    DOUBLE PRECISION, DIMENSION(n1, n2), INTENT(in) :: A
    DOUBLE PRECISION, DIMENSION(n2, n3), INTENT(in) :: B
    DOUBLE PRECISION, DIMENSION(n1, n3), INTENT(out) :: C
    DOUBLE PRECISION, INTENT(out) :: WT
    double precision :: start_time, end_time
    INTEGER :: i, j, k

    C(:,:) = 0.0

    CALL cpu_time(start_time)

    DO j = 1, n3
       DO k = 1, n2
          DO i = 1, n1
             C(i,j) = C(i,j) + A(i,k) * B(k,j)
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_time(end_time)

    WT = end_time - start_time

  END SUBROUTINE jki_mat_mult

  SUBROUTINE f90_mat_mult(n1, n2, n3, A, B, C, WT)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n1, n2, n3
    DOUBLE PRECISION, DIMENSION(n1, n2), INTENT(in) :: A
    DOUBLE PRECISION, DIMENSION(n2, n3), INTENT(in) :: B
    DOUBLE PRECISION, DIMENSION(n1, n3), INTENT(out) :: C
    DOUBLE PRECISION, INTENT(out) :: WT
    double precision :: start_time, end_time

    CALL cpu_time(start_time)

    C = matmul(A, B)

    CALL cpu_time(end_time)

    WT = end_time - start_time

  END SUBROUTINE f90_mat_mult


end program matrix_multiply
