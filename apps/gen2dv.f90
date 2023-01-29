! * gen2dv *
!
! This program generates a velocity field with a given amount of nodes.
! Velocity vertices are to be interpolated by cubic B-splines. The program
! can also add, if required, a random structure and/or a checker board pattern.
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

program gen2dv

  use my_functions
  use my_types
  use m_inout

  implicit none
  type(scalar_field_2d) :: vel
  type(scalar_field_2d) :: cov
  type(scalar_field_2d)	:: ani
  type(gen2dv_conf)     :: conf
  type(list_ind_val)    :: spikes

  real(kind=dbl)::chvp,chvp1,chvp2
  integer ::vusp1,vusp2,vusp1o,vusp2o
  integer i,j

  character(len=strlen)         :: arg


  call get_command_argument(1,arg)
  call read_gen2dv_conf(trim(arg),conf,spikes)

  ! setting up the velocity model and uncertainty grid

  vel%x0=conf%x0;vel%y0=conf%y0
  vel%nx=conf%nx;vel%ny=conf%ny
  vel%dx=conf%dx;vel%dy=conf%dy
  vel%cn=conf%cn
  allocate(vel%val(vel%nx+vel%cn*2,vel%ny+vel%cn*2))
  vel%val=conf%vbg

  ! add checkerboard pattern if required
  if (conf%do_cebo==1) then
     chvp1=conf%cebo_pert
     chvp2=conf%cebo_pert
     chvp=conf%cebo_pert
     vusp1=-1
     vusp2=-1

     do i=1,vel%nx+vel%cn*2
        if (mod(i,conf%cebo_size)==0) then
           chvp1=-chvp1
           if (conf%cebo_spac==1) then
              if(vusp1.EQ.0)THEN
                 if(vusp1o.EQ.-1)THEN
                    vusp1=1
                 else
                    vusp1=-1
                 end if
              else
                 vusp1o=vusp1
                 vusp1=0
              end if
           end if
        end if
        chvp2=chvp1
        vusp2=1
        vusp2o=1

        do j=1,vel%ny+vel%cn*2
           if(MOD(j,conf%cebo_size).EQ.0)THEN
              chvp2=-chvp2
              if(conf%cebo_spac.EQ.1)THEN
                 if(vusp2.EQ.0)THEN
                    if(vusp2o.EQ.-1)THEN
                       vusp2=1
                    else
                       vusp2=-1
                    end if
                 else
                    vusp2o=vusp2
                    vusp2=0
                 end if
              end if
           end if
           chvp=chvp2
           vel%val(i,j)=vel%val(i,j)+vusp1*vusp2*chvp
        end do
     end do
  end if

  ! random perturbations
  if (conf%do_rand==1) then
     do i=1,vel%nx+vel%cn*2
	do j=1,vel%ny+vel%cn*2
           vel%val(i,j)= vel%val(i,j)+gasdev(conf%rseed)*conf%stadev
        end do
     end do
  end if

  ! apply spikes
  do i=1,spikes%n
     vel%val(spikes%ind(i,1),spikes%ind(i,2))=vel%val(spikes%ind(i,1),&
          spikes%ind(i,2))+spikes%val(i)
  end do

  ! write the output to the file
  call write_scalar_field_2d(vel,conf%ofn_velmod)
  call deallocate_scalar_field_2d(vel)


  ! generate uncertainty field (model covariance) if required
  if (conf%do_modcov==1) then
     cov%x0=conf%x0;cov%y0=conf%y0
     cov%nx=conf%nx;cov%ny=conf%ny
     cov%dx=conf%dx;cov%dy=conf%dy
     cov%cn=conf%cn
     allocate(cov%val(cov%nx+cov%cn*2,cov%ny+cov%cn*2))
     cov%val=conf%covbg
     call write_scalar_field_2d(cov,conf%ofn_modcov)
  end if

end program gen2dv

