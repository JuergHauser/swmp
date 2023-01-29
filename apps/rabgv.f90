! * rabgv *
!
! This program resets the velocity values for nodes for which the
! ray path coverage is smaller than a certain threshold.
!
!   This file is part of mps.
!   Copyright (C) 2008 Juerg Hauser
!
!   mps is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   mps is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with mps. If not, see <http://www.gnu.org/licenses/>.
!
!   You may contact the author by:
!       e-mail:  juerg@rses.anu.edu.au
!

program rabgv

  use my_functions
  use my_types
  use m_inout
  implicit none
  type(rabgv_conf):: conf
  type(scalar_field_2d) :: mod
  type(scalar_field_2d) :: rac

  integer i,j
  real(kind=dbl)::dv1,dv2
  character(len=strlen)     :: arg

  call get_command_argument(1,arg)
  call read_rabgv_conf(trim(arg),conf)
  call read_scalar_field_2d(mod,conf%ifn_mod)
  call read_scalar_field_2d(rac,conf%ifn_rac)

  do i=1,mod%nx+mod%cn*2
     do j=1,mod%ny+mod%cn*2
        if (rac%val(i,j)<=conf%mnr) then
           mod%val(i,j)=conf%bgv
        end if

     end do
  end do

  call write_scalar_field_2d(mod,conf%ofn_mod)

end program rabgv

