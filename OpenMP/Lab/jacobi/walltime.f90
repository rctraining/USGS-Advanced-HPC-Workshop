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
! Module to call gettimeofday
!
module walltime

        use, intrinsic :: iso_c_binding, only : c_double

        implicit none

        private
        public  :: get_walltime

        interface
                subroutine get_walltime(wcTime) &
                        bind(c, name="get_walltime")
                        import :: c_double
                        real(c_double), intent(out) :: wcTime
                end subroutine get_walltime
        end interface

end module walltime
