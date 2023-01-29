! * gcray  *
!
! Travel times computation along great circle paths between sources and
! receivers for a constant velocity and a spherical earth
!
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

program gcray

  use my_types
  use my_constants
  use my_functions
  use m_inout

  implicit none

  ! Data strcutures used
  type(receivers)           :: recs
  type(sources)             :: sous
  type(gcray_conf)          :: conf
  type(node),pointer        :: head=>null()
  type(node),pointer        :: tail=>null()

  ! local variables
  integer:: i,j
  real(kind=dbl):: d
  character(len=strlen)     :: arg

  call get_command_argument(1,arg)



  ! read in the configuration and allocate the relevant structures
  call read_gcray_conf(trim(arg),sous,recs,conf)

  ! loop over sources and receivers
  do i=1,sous%nn
     do j=1,recs%nn
        d= great_circle_distance(sous%pos(i,1),recs%pos(j,1),sous%pos(i,2),&
             recs%pos(j,2))
        recs%nar(j)=1
        recs%ariv(j,1)%paic=1
        recs%ariv(j,1)%time=d/conf%vel
     end do
     call write_arrival_times(recs,conf%ofn_traveltimes,i)
  end do

  deallocate(recs%pos,recs%nar,recs%ariv)
  call deallocate_sources(sous)

end program gcray
