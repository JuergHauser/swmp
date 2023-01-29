! * resa2d *
!
! This program reads in a B-spline grid file and dices it up
! for output and plotting using GMT. It is possible to compute
! the gradient of the field
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

program resa2dv

  use my_types
  use my_functions
  use m_inout
  implicit none
  type(scalar_field_2d):: vin,vout,gout
  type(resa2d_conf) :: conf
  integer :: i,j
  real(kind=dbl) :: x,y,f,fx,fy
  character(len=strlen)     :: arg

  call get_command_argument(1,arg)
  call read_resa2d_conf(trim(arg),conf)

  ! load the velocity fiel conf%igrid

  open(unit=input,file=conf%igrid,status='old')
  read(input,*) vin%x0,vin%y0
  read(input,*) vin%nx,vin%ny
  read(input,*) vin%dx,vin%dy
  read(input,*) vin%cn
  allocate(vin%val(vin%nx+2*vin%cn,vin%ny+2*vin%cn))
  do i=1,vin%nx+2*vin%cn
     do j=1,vin%ny+2*vin%cn
        read(input,*) vin%val(i,j)
     end do
  end do

  ! setup the sacalar field 2d grid for the updiced velocity field and it's
  ! gradient
  vout%x0=vin%x0;vout%y0=vin%y0
  vout%nx=conf%nrx*(vin%nx-1)+1;vout%ny=conf%nry*(vin%ny-1)+1
  vout%dx=((vin%nx-1)*vin%dx)/(vout%nx-1)
  vout%dy=((vin%ny-1)*vin%dy)/(vout%ny-1)
  vout%cn=vin%cn
  allocate(vout%val(vout%nx+vout%cn*2,vout%ny+vout%cn*2))

  gout%x0=vin%x0;gout%y0=vin%y0
  gout%nx=conf%nrx*(vin%nx-1)+1;gout%ny=conf%nry*(vin%ny-1)+1
  gout%dx=((vin%nx-1)*vin%dx)/(gout%nx-1)
  gout%dy=((vin%ny-1)*vin%dy)/(gout%ny-1)
  gout%cn=vin%cn
  allocate(gout%val(gout%nx+gout%cn*2,gout%ny+gout%cn*2))

  ! compute the velocit field and its gradient for the nodes of the
  ! updiced grid using a cubic bspline interpolation
  do i=1+vout%cn,vout%nx+vout%cn
     do j=1+vout%cn,vout%ny+vout%cn
        x=(i-1-vout%cn)*vout%dx+vout%x0
        y=(j-1-vout%cn)*vout%dy+vout%y0

        call cubic_bspline_interp_val_1der(vin,x,y,f,fx,fy)
        vout%val(i,j)=f
        gout%val(i,j)=sqrt(fx**2+fy**2)
     end do
  end do

  ! write velocity field to files
  if (conf%sw==1.or.conf%sw==3) then
     call write_gmt_scalar_field_2d(vout,conf%ovef,conf%gmtvef)
  end if

  if (conf%sw==2.or.conf%sw==3) then
     call write_gmt_scalar_field_2d(gout,conf%ogaf,conf%gmtgaf)
  end if

  call deallocate_scalar_field_2d(vin)
  call deallocate_scalar_field_2d(vout)
  call deallocate_scalar_field_2d(gout)

end program resa2dv
