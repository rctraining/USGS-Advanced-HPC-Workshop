!
! Copyright (C) 2016  Timothy Brown
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!
! Module for error message reporting.
!

module err

        implicit none

        private

        integer, parameter :: EX_SOFTWARE = 70

        public   :: err_msg,  &
                    err_abort

contains

!
! Write an error message to stderr.
!
subroutine err_msg(message)

        use iso_fortran_env, only : error_unit

        implicit none

        character(len=*), intent(in) :: message

        write(error_unit,*) trim(message)

end subroutine err_msg

!
! Write the error message using err_msg() and then stop the program.
!
subroutine err_abort(message, ierr)

        implicit none

        character(len=*),  intent(in) :: message
        integer, optional, intent(in) :: ierr

        integer :: ecode

        call err_msg(message)

        if (present(ierr)) then
                ecode = ierr
        else
                ecode = EX_SOFTWARE
        end if

        ! ifort does not like having a variable given to stop.
        stop 70

end subroutine err_abort

end module err
