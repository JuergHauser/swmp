!
!   This file is part of swmp.
!   Copyright (C) 2023 CSIRO
!
!   swmp is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   swmp is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with swmp. If not, see <http://www.gnu.org/licenses/>.
!
!   You may contact the author by:
!       e-mail:  juerg.hauser@csiro.au
!

! * m_error *
!
! Error handling module for C/Python interface
!

module m_error
  use, intrinsic :: iso_c_binding
  implicit none

  ! Module-level error state
  integer :: error_code = 0
  character(256) :: error_message = ""

  ! Error code constants
  integer, parameter :: SWMP_SUCCESS = 0
  integer, parameter :: SWMP_ERROR_FILE = -1
  integer, parameter :: SWMP_ERROR_MEMORY = -2
  integer, parameter :: SWMP_ERROR_INVALID_PARAM = -3

contains

  !---------------------------------------------------------------------------!

  subroutine set_error(code, msg)
    ! Set error state (internal use)

    implicit none
    integer, intent(in) :: code
    character(*), intent(in) :: msg

    error_code = code
    error_message = trim(msg)

  end subroutine set_error

  !---------------------------------------------------------------------------!

  function get_last_error(code, msg_ptr, msg_len) result(status) bind(c, name="get_last_error")
    ! Get last error from C/Python
    !
    ! Parameters:
    !   code (out) - error code
    !   msg_ptr (in) - pointer to message buffer
    !   msg_len (in) - buffer length
    !
    ! Returns:
    !   status - always 0 (success)

    implicit none
    integer(c_int), intent(out) :: code
    type(c_ptr), value :: msg_ptr
    integer(c_int), value :: msg_len
    integer(c_int) :: status

    character(len=msg_len, kind=c_char), pointer :: msg_str

    code = error_code
    call c_f_pointer(msg_ptr, msg_str)
    msg_str = trim(error_message)
    status = SWMP_SUCCESS

  end function get_last_error

  !---------------------------------------------------------------------------!

  function clear_error() result(status) bind(c, name="clear_error")
    ! Clear error state from C/Python
    !
    ! Returns:
    !   status - always 0 (success)

    implicit none
    integer(c_int) :: status

    error_code = 0
    error_message = ""
    status = SWMP_SUCCESS

  end function clear_error

  !---------------------------------------------------------------------------!

end module m_error
