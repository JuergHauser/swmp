! * inv *
!
! Non linear inversion part of mps using the pseudo inverse or the
! subspace method.
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

program inv

  use my_types
  use my_functions
  use m_inout
  use m_linalg
  use m_inverse
  implicit none

  type(inv_conf)	   :: conf

  type(observed_arrivals)  :: obs
  type(predicted_arrivals) :: pred

  type(scalar_field_2d)    :: vcur
  type(scalar_field_2d)    :: vpert
  type(scalar_field_2d)    :: rapco  ! ray path coverage
  type(list_index)         :: usarr  ! list of rays used

  ! types for the inversion
  type(vector)     ::  dm      ! model perturbation
  type(vector)     ::  m       ! current model
  type(vector)     ::  m0      ! starting model
  type(matrix)     ::  a       ! projection matrix
  type(crs_matrix) ::  g       ! frechet matrix
  type(crs_matrix) ::  gt      ! transpose of the frechet matrix
  type(crs_matrix) ::  cdi     ! inverse of the data covariance matrix
  type(crs_matrix) ::  cm      ! model covariance matrix
  type(crs_matrix) ::  cmi     ! inverse of the model covariance matrix
  type(crs_matrix) ::  dtd     ! smoothing operator
  type(vector)     ::  dobs    ! observation
  type(vector)     ::  dpred   ! model prediction
  type(vector)     ::  gamma   ! gamma
  real(kind=dbl)   ::  epsilon ! damping factor
  real(kind=dbl)   ::  eta     ! smoothing factor


  character(len=strlen)         :: arg

  ! read in the configurations
  call get_command_argument(1,arg)
  call read_inv_conf(trim(arg),conf,epsilon,eta)

  ! read predicted and observ ed traveltimes
  call read_predicted_arrivals(pred,conf%ifn_predtimes)
  call read_observed_arrivals(obs,conf%ifn_obstimes)

  ! generate dobs and dpred
  call initialise_data(obs,pred,dobs,cdi,dpred,usarr)

  ! create the frechet matrix the observations and the model prediction vector
  call read_frechet_matrix(g,usarr,conf%ifn_frechet_mat,conf%ifn_frechet_hdr,conf%ifn_frechet_rai)

  ! get the velocity model and the model covariance
  call initialise_model(m0,vcur,m,cm,cmi,conf)

  select case (conf%invmod)
  case(1)
     call compute_smoothing_operator(vcur,m,dtd,conf)
     call compute_projection_matrix(g,gt,dobs,cdi,dpred,m0,cm,cmi,m,dtd,&
          epsilon,eta,gamma,a,conf)
     call compute_model_perturbation(g,gt,cdi,cmi,dtd,gamma,a,&
          epsilon,eta,dm,conf)
  case(2)
     call use_pseudo_inverse(g,dobs,dpred,dm,conf)
  end select

  ! compute raypath coverage
  call compute_raypath_coverage(g,vcur,rapco)

  ! update the velocity mode
  call update_velocity_model(vcur,vpert,dm)

  ! write the new velocity model
  call write_scalar_field_2d(vpert,conf%ofn_vpert)
  call write_scalar_field_2d(vcur,conf%ofn_vupd)

  ! write raypath coverage
  call write_scalar_field_2d(rapco,conf%ofn_rapco)

! cleanup
call deallocate_vector(dm)
call deallocate_vector(m)
call deallocate_vector(m0)
call deallocate_matrix(a)
call deallocate_crs_matrix(g)
call deallocate_crs_matrix(gt)
call deallocate_crs_matrix(cdi)
call deallocate_crs_matrix(cm)
call deallocate_crs_matrix(cmi)
call deallocate_crs_matrix(dtd)
call deallocate_vector(dobs)
call deallocate_vector(dpred)
call deallocate_vector(gamma)

call deallocate_scalar_field_2d(vcur)
call deallocate_scalar_field_2d(vpert)
call deallocate_scalar_field_2d(rapco)

call deallocate_list_index(usarr)

call deallocate_observed_arrivals(obs)
call deallocate_predicted_arrivals(pred)

end program inv
