PROGRAM MAIN
    USE ISO_FORTRAN_ENV, ONLY : output_unit
    USE MPI
    IMPLICIT NONE
    INTEGER :: ncpu, my_rank, ierr

    CALL INIT_COMM()

    CALL SWAP_DOUBLE()

    CALL SWAP_INTEGER()

    CALL SWAP_MULTI_DOUBLE()

    CALL MPI_FINALIZE(ierr)

CONTAINS

    SUBROUTINE SWAP_DOUBLE()
        IMPLICIT NONE
        REAL*8 :: my_val, send_val
        INTEGER :: mtag1 = 3, mtag2 = 4, ierr
        INTEGER :: source, dest, ireq(2)
        INTEGER :: mstat(MPI_STATUS_SIZE, 2)
        my_val = my_rank
        !Swap the values of my_val between rank 0 and 1
        send_val = my_val


        IF (my_rank .eq. 0) THEN 

            source = 1  ! receive from rank 1
            ! The "1" in the calls to receive and send indicates that we wish to transmit
            ! a single value of type MPI_DOUBLE_PRECISION (Real*8)
            WRITE(output_unit,*)' '
            WRITE(output_unit,*)'     Sending/Receiving one double-precision (8-byte) number...'
            WRITE(output_unit,'(a,F4.1)')'Rank 0;  pre-receive; my_val = ', my_val
            Call MPI_irecv(my_val, 1, MPI_DOUBLE_PRECISION, source, & 
                mtag1, MPI_COMM_WORLD, ireq(1),ierr)


            !Once we have received, we send to rank 1
            dest = 1
            Call MPI_isend(send_val, 1, MPI_DOUBLE_PRECISION, dest, & 
                mtag2, MPI_COMM_WORLD, ireq(2),ierr)
        ENDIF

        IF (my_rank .eq. 1) THEN ! highest rank does not send
           
            ! Mimic the logic used for rank 0, but send first and swap the tags
            dest = 0
            Call MPI_isend(send_val, 1, MPI_DOUBLE_PRECISION, dest, & 
                mtag1, MPI_COMM_WORLD, ireq(1),ierr)

            source = 0 
            WRITE(output_unit,'(a,F4.1)')'Rank 1;  pre-receive; my_val = ', my_val
            Call MPI_irecv(my_val, 1, MPI_DOUBLE_PRECISION, source, & 
                mtag2, MPI_COMM_WORLD, ireq(2),ierr)


        ENDIF
        CALL MPI_Waitall(2, ireq, mstat, ierr)
        WRITE(output_unit,'(a,i0,a,F4.1)')'Rank ',my_rank,'; post-receive; my_val = ', my_val
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  ! Pause execution to give I/O time to catch up
    END SUBROUTINE SWAP_DOUBLE

    SUBROUTINE SWAP_INTEGER()
        IMPLICIT NONE
        INTEGER :: my_val, send_val
        INTEGER :: mtag1 = 3, mtag2 = 4, ierr
        INTEGER :: source, dest, ireq(2)
        INTEGER :: mstat(MPI_STATUS_SIZE, 2)
        my_val = my_rank
        !Swap the values of my_val between rank 0 and 1
        send_val = my_val


        IF (my_rank .eq. 0) THEN 

            source = 1  ! receive from rank 1
            ! The "1" in the calls to receive and send indicates that we wish to transmit
            ! a single value of type MPI_DOUBLE_PRECISION (Real*8)
            WRITE(output_unit,*)' '
            WRITE(output_unit,*)'     Sending/Receiving one integer (4-byte) number...'
            WRITE(output_unit,'(a,i0)')'Rank 0;  pre-receive; my_val = ', my_val
            Call MPI_irecv(my_val, 1, MPI_INTEGER, source, & 
                mtag1, MPI_COMM_WORLD, ireq(1),ierr)


            !Once we have received, we send to rank 1
            dest = 1
            Call MPI_isend(send_val, 1, MPI_INTEGER, dest, & 
                mtag2, MPI_COMM_WORLD, ireq(2),ierr)
        ENDIF

        IF (my_rank .eq. 1) THEN ! highest rank does not send
           
            ! Mimic the logic used for rank 0, but send first and swap the tags
            dest = 0
            Call MPI_isend(send_val, 1, MPI_INTEGER, dest, & 
                mtag1, MPI_COMM_WORLD, ireq(1),ierr)

            source = 0 
            WRITE(output_unit,'(a,i0)')'Rank 1;  pre-receive; my_val = ', my_val
            Call MPI_irecv(my_val, 1, MPI_INTEGER, source, & 
                mtag2, MPI_COMM_WORLD, ireq(2),ierr)


        ENDIF
        CALL MPI_Waitall(2, ireq, mstat, ierr)
        WRITE(output_unit,'(a,i0,a,i0)')'Rank ',my_rank,'; post-receive; my_val = ', my_val
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  ! Pause execution to give I/O time to catch up
    END SUBROUTINE SWAP_INTEGER

    SUBROUTINE SWAP_MULTI_DOUBLE()
        IMPLICIT NONE
        REAL*8 :: my_vals(1:4), send_buffer(1:2)
        INTEGER :: mtag1 = 3, mtag2 = 4, ierr
        INTEGER :: source, dest, ireq(2)
        INTEGER :: mstat(MPI_STATUS_SIZE, 2)
        my_vals(:) = my_rank
        !Swap the last two elements of my_vals between ranks 0 and 1
        send_buffer(1:2) = my_vals(3:4)


        IF (my_rank .eq. 0) THEN 

            ! Receive from rank 1
            ! While we could receive into a separate receive buffer,
            ! we can also receive directly into the vals array.
            ! Note that this wouldn't work reliably if we did not first
            ! copy the elements we want to send into the send_buffer.

            !In this case, we look for a message containing two double-precision
            ! elements and store them in memory beginning at my_vals(3)
            !    -- had we not specified my_vals(3), and just my_vals, say, 
            !       data received would be copied into memory beginning at my_vals(1).
            source = 1  ! receive from rank 1
            WRITE(output_unit,*)' '
            WRITE(output_unit,*)'     Sending/Receiving multiple double-precision numbers...'
            WRITE(output_unit,'(a,4F4.1)')'Rank 0;  pre-receive; my_vals = ', my_vals
            Call MPI_irecv(my_vals(3), 2, MPI_DOUBLE_PRECISION, source, & 
                mtag1, MPI_COMM_WORLD, ireq(1),ierr)

            !Once we have received, we send to rank 1
            dest = 1
            Call MPI_isend(send_buffer, 2, MPI_DOUBLE_PRECISION, dest, & 
                mtag2, MPI_COMM_WORLD, ireq(2),ierr)
        ENDIF
        IF (my_rank .eq. 1) THEN ! highest rank does not send
           
            ! Mimic the logic used for rank 0, but send first and swap the tags.
            ! Q:  What happens if we try to receive first?


            dest = 0
            Call MPI_isend(send_buffer, 2, MPI_DOUBLE_PRECISION, dest, & 
                mtag1, MPI_COMM_WORLD, ireq(1),ierr)

            source = 0 
            WRITE(output_unit,'(a,4F4.1)')'Rank 1;  pre-receive; my_vals = ', my_vals
            Call MPI_irecv(my_vals(3), 2, MPI_DOUBLE_PRECISION, source, & 
                mtag2, MPI_COMM_WORLD, ireq(2),ierr)


        ENDIF
        CALL MPI_Waitall(2, ireq, mstat, ierr)
        WRITE(output_unit,'(a,i0,a,4F4.1)')'Rank ',my_rank,'; post-receive; my_val = ', my_vals
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  ! Pause execution to give I/O time to catch up
    END SUBROUTINE SWAP_MULTI_DOUBLE

    SUBROUTINE INIT_COMM()
        IMPLICIT NONE
        INTEGER :: ierr
        !Initialize MPI  (ierr is an error flag with nonzero value if there's a problem)
        Call MPI_INIT( ierr )
        IF (ierr .ne. 0) STOP 'Error intializing MPI'

        !Find number of MPI processes executing this program (ncpu)
        Call MPI_Comm_size(MPI_COMM_WORLD, ncpu,ierr)
        IF (ierr .ne. 0) STOP 'Error finding ncpu'

        !Find this process's rank (my_rank; unique numeric identifier for each MPI process)
        Call MPI_Comm_rank(MPI_COMM_WORLD, my_rank,ierr)
        IF (ierr .ne. 0) STOP 'Error finding my_rank'
        WRITE(output_unit,'(a,i0,a,i0,a)')'Rank ',my_rank,' of ', ncpu, ' is active.'

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  ! Pause until all processes have initiated MPI

    END SUBROUTINE INIT_COMM
END PROGRAM 
