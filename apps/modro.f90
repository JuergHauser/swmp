! * per2dv *
!
! This program calcualtes an estimate of model roughness and variance
! by dicing up a cubic B-spline grid
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

program modro

  use my_functions
  use my_types
  use my_constants
  use m_inout

  implicit none
  type(scalar_field_2d)      :: vcur
  type(scalar_field_2d)      :: vref
  type(scalar_field_2d)      :: vups
  type(modro_conf)           :: conf
  real(kind=dbl)             :: mrough,mvar
  integer :: i,j
  real(kind=dbl)::x,y,fref,fcur
  real(kind=dbl)::ri,risti,dp,dt
  character(len=strlen)     :: arg

  call get_command_argument(1,arg)
  call read_modro_conf(trim(arg),conf)

  call read_scalar_field_2d(vcur,conf%ifn_vcur)
  call read_scalar_field_2d(vref,conf%ifn_vref)

  vups%x0=vcur%x0;vups%y0=vcur%y0
  vups%nx=conf%nrx*(vcur%nx-1)+1;vups%ny=conf%nry*(vcur%ny-1)+1
  vups%dx=((vcur%nx-1)*vcur%dx)/(vups%nx-1)
  vups%dy=((vcur%ny-1)*vcur%dy)/(vups%ny-1)
  vups%cn=vcur%cn
  allocate(vups%val(vups%nx+vups%cn*2,vups%ny+vups%cn*2))

mvar=0.0
  do i=1+vups%cn,vups%nx+vups%cn
     do j=1+vups%cn,vups%ny+vups%cn
        x=(i-1-vups%cn)*vups%dx+vups%x0
        y=(j-1-vups%cn)*vups%dy+vups%y0
        call cubic_bspline_interp_val(vcur,x,y,fref)
        call cubic_bspline_interp_val(vref,x,y,fcur)
        vups%val(i,j)=fcur-fref
        mvar=mvar+(fcur-fref)**2
     end do
  end do

  mvar=mvar/(vups%nx*vups%ny-1.0)

mrough=0.0
do i=1+vups%cn+1,vups%nx+vups%cn-2
   do j=1+vups%cn+1,vups%ny+vups%cn-2
      ri=conf%rearth
      risti=ri*sin((90.0-vcur%y0)*pi/180.0+(j-1)*vups%dy)
      dp=vups%val(i,j+1)-2.0*vups%val(i,j)+vups%val(i,j-1)
      dp=dp/((risti*vups%dx)**2)
      dt=vups%val(i+1,j)-2.0*vups%val(i,j)+vups%val(i-1,j)
      dt=dt/((ri*vups%dy)**2)
      mrough=mrough+abs(dp)+abs(dt)
   end do
end do
mrough=mrough/((vups%nx-2)*(vups%ny-2))

WRITE(6,*)'Model variance in (km/s)**2 is ',mvar
WRITE(6,*)'Model roughness in (kms)**(-1) is ',mrough



end program modro

