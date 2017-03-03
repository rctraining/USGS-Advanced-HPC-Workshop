program simd_f
        implicit none

        integer, parameter            :: n = 64000000

        integer                       :: ierr
        double precision, allocatable :: a(:)
        double precision, allocatable :: b(:)
        double precision, allocatable :: c(:)
        !dir$ attributes align: __BIGGEST_ALIGNMENT__ :: a
        !dir$ attributes align: __BIGGEST_ALIGNMENT__ :: b
        !dir$ attributes align: __BIGGEST_ALIGNMENT__ :: c

        ierr = 0
        allocate(a(n), stat=ierr)
        if (ierr .ne. 0) then
                write(0,*) 'Unable to allocate: ', n
                stop
        end if
        allocate(b(n), stat=ierr)
        if (ierr .ne. 0) then
                write(0,*) 'Unable to allocate: ', n
                stop
        end if
        allocate(c(n), stat=ierr)
        if (ierr .ne. 0) then
                write(0,*) 'Unable to allocate: ', n
                stop
        end if

        a = 1
        b = 2
        c = 0

        call plus(a, b, c)

!        write(6, *) 'a: ', a
!        write(6, *) 'b: ', b
!        write(6, *) 'c: ', c

        if (allocated(a)) then
                deallocate(a)
        end if
        if (allocated(b)) then
                deallocate(b)
        end if
        if (allocated(c)) then
                deallocate(c)
        end if
contains

subroutine plus(a, b, c)
        implicit none
        double precision, intent(in   ) :: a(:)
        double precision, intent(in   ) :: b(:)
        double precision, intent(  out) :: c(:)

        integer                         :: i
        integer                         :: n

        !dir$ assume_aligned a: __BIGGEST_ALIGNMENT__
        !dir$ assume_aligned b: __BIGGEST_ALIGNMENT__
        !dir$ assume_aligned c: __BIGGEST_ALIGNMENT__

        n = ubound(a, 1)

!$OMP parallel
!$OMP do simd &
!$OMP aligned(a,b,c:__BIGGEST_ALIGNMENT__)
        do i=1,n
                c(i) = a(i) + b(i)
        end do
!$OMP end do simd
!$OMP end parallel
end subroutine plus

end program simd_f
