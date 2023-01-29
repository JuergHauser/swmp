! * tresid *
!
! This program analyses the travel time residuals
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
!   along with mps. if not, see <http://www.gnu.org/licenses/>.
!
!   You may contact the author by:
!       e-mail:  juerg@rses.anu.edu.au
!

program tresid

  use my_functions
  use my_types
  use m_inout
  use m_inverse
  implicit none

  type(observed_arrivals)  :: obs
  type(predicted_arrivals) :: pred
  type(traveltime_residuals):: resid
  type(tresid_conf):: conf
  character(len=strlen)     :: arg

  call get_command_argument(1,arg)
  call read_tresid_conf(trim(arg),conf)

  ! read predicted and observed traveltimes
  call read_predicted_arrivals(pred,conf%ifn_predtimes)
  call read_observed_arrivals(obs,conf%ifn_obstimes)

  ! compute residuals
  call compute_residuals(obs,pred,resid)

  ! write results
  call write_tresid_summary(resid,conf)
  call write_residuals(resid,conf)

end program tresid

