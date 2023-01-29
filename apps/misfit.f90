! * misfit *
!
! This program comptes the difference between two scalar fields
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

program misfit

  use my_types
  use my_functions
  use m_inout
  implicit none
  type(scalar_field_2d):: vrefi,vtruei,vinvi
   type(scalar_field_2d):: vref,vtrue,vinv

  type(scalar_field_2d):: atrue,ainv
  type(scalar_field_2d):: dabs,drel
  type(scalar_field_2d):: adrel

  real(kind=dbl):: rms_abs,cor_abs(3),cor_ano(3)
  real(kind=dbl),allocatable::v1(:),v2(:)
  integer :: i,j,k,n
  character(len=strlen)     :: arg

  type(misfit_conf) :: conf

   call get_command_argument(1,arg)

  call read_misfit_conf(trim(arg),conf)

  call read_scalar_field_2d(vinvi,conf%ifn_vinv)
  call read_scalar_field_2d(vtruei,conf%ifn_vtrue)
  call read_scalar_field_2d(vrefi,conf%ifn_vref)

  ! updice the velocity grids
  call upsample_bspline(vinvi,vinv,conf%nrx,conf%nry)
  call upsample_bspline(vtruei,vtrue,conf%nrx,conf%nry)
  call upsample_bspline(vrefi,vref,conf%nrx,conf%nry)

  ! compute absolute difference grid
  dabs%x0=vref%x0;dabs%y0=vref%y0
  dabs%nx=vref%nx;dabs%ny=vref%ny
  dabs%dx=vref%dx;dabs%dy=vref%dy
  dabs%cn=vref%cn
  allocate(dabs%val(dabs%nx+dabs%cn*2,dabs%ny+dabs%cn*2))
  do i=1,dabs%nx+dabs%cn*2
     do j=1,dabs%ny+dabs%cn*2
        dabs%val(i,j)=vtrue%val(i,j)-vinv%val(i,j)
     end do
  end do

  drel%x0=vref%x0;drel%y0=vref%y0
  drel%nx=vref%nx;drel%ny=vref%ny
  drel%dx=vref%dx;drel%dy=vref%dy
  drel%cn=vref%cn
  allocate(drel%val(drel%nx+drel%cn*2,drel%ny+drel%cn*2))
  ! compute relative difference grid
  do i=1,drel%nx+drel%cn
     do j=1,drel%ny+drel%cn
        drel%val(i,j)=(vtrue%val(i,j)-vinv%val(i,j))/vtrue%val(i,j)*100.0_dbl
     end do
  end do

 ! compute the anomalies

  atrue%x0=vref%x0;atrue%y0=vref%y0
  atrue%nx=vref%nx;atrue%ny=vref%ny
  atrue%dx=vref%dx;atrue%dy=vref%dy
  atrue%cn=vref%cn
  allocate(atrue%val(atrue%nx+atrue%cn*2,atrue%ny+atrue%cn*2))
  do i=1,atrue%nx+atrue%cn*2
 do j=1,atrue%ny+atrue%cn*2
 atrue%val(i,j)=vtrue%val(i,j)-vref%val(i,j)
 end do
 end do

  ainv%x0=vref%x0;ainv%y0=vref%y0
  ainv%nx=vref%nx;ainv%ny=vref%ny
  ainv%dx=vref%dx;ainv%dy=vref%dy
  ainv%cn=vref%cn
  allocate(ainv%val(ainv%nx+ainv%cn*2,ainv%ny+ainv%cn*2))
  do i=1,ainv%nx+ainv%cn*2
 do j=1,ainv%ny+ainv%cn*2
 ainv%val(i,j)=vinv%val(i,j)-vref%val(i,j)
 end do
 end do


 ! compute the relative difference between the anomalies
  adrel%x0=vref%x0;adrel%y0=vref%y0
  adrel%nx=vref%nx;adrel%ny=vref%ny
  adrel%dx=vref%dx;adrel%dy=vref%dy
  adrel%cn=vref%cn
  allocate(adrel%val(adrel%nx+adrel%cn*2,adrel%ny+adrel%cn*2))
  do i=1,adrel%nx+adrel%cn*2
     do j=1,adrel%ny+adrel%cn*2
        adrel%val(i,j)=(atrue%val(i,j)-ainv%val(i,j))/atrue%val(i,j)*100.0_dbl
     end do
  end do


! compute the two rms errors

   rms_abs=0.0

  do i=1,dabs%nx+dabs%cn*2
     do j=1,dabs%ny+dabs%cn*2
        rms_abs=rms_abs+(dabs%val(i,j))**2
     end do
  end do

  rms_abs=sqrt(rms_abs/((dabs%nx+dabs%cn*2)*(dabs%ny+dabs%cn*2)))


 ! compute the correlation coefficient
  n=((vref%nx+vref%cn*2)*(vref%ny+vref%cn*2))
  allocate (v1(n))
  allocate(v2(n))

  k=0

  do i=1,vref%nx+vref%cn*2
     do j=1,vref%ny+vref%cn*2
        k=k+1
        v1(k)=vinv%val(i,j)
        v2(k)=vtrue%val(i,j)
     end do
  end do

  call pearsn(v1,v2,cor_abs(1),cor_abs(2),cor_abs(3))


  k=0

  do i=1,vref%nx+vref%cn*2
     do j=1,vref%ny+vref%cn*2
        k=k+1
        v1(k)=ainv%val(i,j)
        v2(k)=atrue%val(i,j)
     end do
  end do

  call pearsn(v1,v2,cor_ano(1),cor_ano(2),cor_ano(3))

 !!  write(*,'(A,f8.5)') 'confidence level        ', 1.0_dbl-prob  prob=cor_ano(2)

!!write(*,*) cor_ano(2)

call write_scalar_field_2d(dabs,conf%ofn_dabs)
call write_scalar_field_2d(drel,conf%ofn_drel)
call write_scalar_field_2d(adrel,conf%ofn_adrel)
call write_misfit_summary(rms_abs,cor_abs,cor_ano,conf)


end program misfit
