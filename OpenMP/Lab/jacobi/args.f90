!
! Module for parsing command line args
!

module args

        use kinds,  only : r_dp
        use err,    only : err_msg

        implicit none

        private
        public   :: args_parse, &
                    args_usage

contains

subroutine args_parse(i, j, k, niter,                        &
                      north, south, east, west, top, bottom, &
                      ierr, output)

        implicit none

        integer,                    intent(out)   :: i
        integer,                    intent(out)   :: j
        integer,                    intent(out)   :: k
        integer,                    intent(out)   :: niter
        real(kind=r_dp),            intent(out)   :: north
        real(kind=r_dp),            intent(out)   :: south
        real(kind=r_dp),            intent(out)   :: east
        real(kind=r_dp),            intent(out)   :: west
        real(kind=r_dp),            intent(out)   :: top
        real(kind=r_dp),            intent(out)   :: bottom
        integer,                    intent(out)   :: ierr
        character(len=*), optional, intent(out)   :: output

        integer :: n                    ! Number of args
        integer :: x, y                 ! Argument iterators
        integer :: iunit                ! Namelist file unit
        character(len=1024) :: key      ! Argument key
        character(len=1024) :: val      ! Argument value
        character(len=1024) :: nlfile   ! Namelist file name

        namelist /domain/ i, j, k
        namelist /jacobi/ niter, north, south, east, west, &
                          top, bottom

        ierr = 0
        nlfile = ''
        output = ''
        y = 0
        n = command_argument_count()
        do x=1,n,2
                call get_command_argument(x, key)
                call get_command_argument(x+1, val)
                select case(key)
                        case('-imax')
                                read(val, '(I8)') i
                                y = y + 1
                        case('-jmax')
                                read(val, '(I8)') j
                                y = y + 1
                        case('-kmax')
                                read(val, '(I8)') k
                                y = y + 1
                        case('-niter')
                                read(val, '(I8)') niter
                                y = y + 1
                        case('-east')
                                read(val, '(F8.0)') east
                                y = y + 1
                        case('-west')
                                read(val, '(F8.0)') west
                                y = y + 1
                        case('-north')
                                read(val, '(F8.0)') north
                                y = y + 1
                        case('-south')
                                read(val, '(F8.0)') south
                                y = y + 1
                        case('-top')
                                read(val, '(F8.0)') top
                                y = y + 1
                        case('-bottom')
                                read(val, '(F8.0)') bottom
                                y = y + 1
                        case('-f')
                                nlfile = val
                        case('-o')
                                output = val
                        case default
                                call err_msg('Unrecognized option: ' &
                                             // trim(key))
                end select
        end do

        if (len(trim(nlfile)) > 0) then
                iunit = 7
                open(unit=iunit, file=nlfile, action='read')
                read(iunit, nml=domain)
                read(iunit, nml=jacobi)
                close(iunit)
                y = 10
        end if

        if ( y /= 10) then
                call err_msg('Not enough command line arguments.')
                call args_usage
                ierr = 1
        end if

end subroutine args_parse

subroutine args_usage()

        implicit none

        character(len=1024) :: prog_name

        call get_command_argument(0, prog_name)

        write(0,*) ''
        write(0,*) 'usage: ', trim(prog_name),  &
                   ' -imax X -jmax Y -kmax Z \ '
        write(0,*) '       -north A -south B -east C -west D -top E ', &
                   '-bottom F -niter N \'
        write(0,*) '       [-f namelist] [-o output]'
        write(0,*) ''
        write(0,*) '  -imax   X    Size of domain in the i direction.'
        write(0,*) '  -jmax   Y    Size of domain in the j direction.'
        write(0,*) '  -kmax   Z    Size of domain in the k direction.'
        write(0,*) '  -north  A    Initial condition on the north face.'
        write(0,*) '  -south  B    Initial condition on the south face.'
        write(0,*) '  -east   C    Initial condition on the east face.'
        write(0,*) '  -west   D    Initial condition on the west face.'
        write(0,*) '  -top    E    Initial condition on the top face.'
        write(0,*) '  -bottom F    Initial condition on the bottom face.'
        write(0,*) '  -niter  N    Number of iterations to make.'
        write(0,*) '  -f namelist  Namelist file.'
        write(0,*) '  -o output    HDF5 to write iterations to.'

end subroutine args_usage

end module args
