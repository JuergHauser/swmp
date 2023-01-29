! * creobs *
!
! This program adds random noise and data covariance information to a set of 
! travel times.
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

program creobs

  use my_types
  use my_functions
  use m_inout
  implicit none

  type(creobs_conf)       :: conf
  type(predicted_arrivals) :: pred
  type(observed_arrivals)  :: obs
  integer :: i,j

  ! load necessary files
  call read_creobs_conf(conf)
  call read_predicted_arrivals(pred,conf%ifn_traveltimes)

  ! copy pred into obs
  obs%n=pred%n
  allocate(obs%sou(obs%n))
  allocate(obs%rec(obs%n))
  allocate(obs%arn(obs%n))
  allocate(obs%art(obs%n))
  allocate(obs%azi(obs%n))

  allocate(obs%unc(obs%n))

  obs%sou=pred%sou
  obs%rec=pred%rec
  obs%arn=pred%arn
  obs%art=pred%art
  obs%azi=pred%azi

  ! add random nois to the arrivla times

  if (conf%do_ran==1) then
     do i=1,obs%n
	obs%art(i)=pred%art(i)+gasdev(conf%rseed)*conf%stdev
     end do
  end if

  obs%unc=conf%unc

  call write_observed_arrivals(obs,conf%ofn_traveltimes)

call deallocate_predicted_arrivals(pred)
call deallocate_observed_arrivals(obs)

end program creobs
