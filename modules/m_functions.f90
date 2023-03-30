!
!   This file is part of swmp.
!   Copyright (C) 2023 CSIRO
!
!   swmp is free software: you can redistribute it and/or modify
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
!       e-mail:  juerg.hauser@csiro.au
!
! m_functions
!
! This module contains a set of basic functions and subroutines used
! through out the programs and modules belonging to the mpstomo package
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

module my_functions


use my_types
  implicit none

contains

  !*****************************************************************************!
  ! procedures used for standardised error and warning handling                 !
  !*****************************************************************************!

  subroutine nrerror(proc,err)
    implicit none
    character(len=*), intent(in) :: proc,err
    write (*,*) '*************************************************************'
    write (*,*) 'procedure: ',proc
    write (*,*) err
    write (*,*) 'program terminated due to numerical error'
    write (*,*) '*************************************************************'
    stop
  end subroutine nrerror

  !----------------------------------------------------------------------------!

  subroutine nrwarn(proc,warn)
    implicit none
    character(len=*), intent(in) :: proc,warn
    write (*,*) '*************************************************************'
    write (*,*) 'procedure: ',proc
    write (*,*) warn
    write(*,*) 'program halted due to numerical problem (enter to continue)'
    write (*,*) '*************************************************************'
    read(*,*)
  end subroutine nrwarn

  !-----------------------------------------------------------------------------!

  subroutine generror(proc,err)
    implicit none
    character(len=*), intent(in) :: proc,err
    write (*,*) '*************************************************************'
    write (*,*) 'procedure: ',proc
    write (*,*) err
    write (*,*) 'program terminated due to general error'
    write (*,*) '*************************************************************'
    stop
  end subroutine generror

  !----------------------------------------------------------------------------!

  subroutine genmsg(proc,msg1,msg2,n)
    implicit none
    integer :: n
    character(len=*), intent(in) :: proc,msg1,msg2
    write (*,*) '*************************************************************'
    write (*,*) 'procedure: ',proc
    write (*,*) msg1
    if (n==2) then
       write (*,*) msg2
    end if
    write (*,*) '*************************************************************'
  end subroutine genmsg

  !-----------------------------------------------------------------------------!

  subroutine heap_sort_list_index(ar)
    ! Heap sort algorithm for a list_index.
    ! based on numerical recepies in fortran
    use my_types
    implicit none

    type(list_index):: ar
    integer :: i

    do i=int(ar%n/2),1,-1
       call shift_down(ar,i,ar%n)
    end do
    do i=int(ar%n),2,-1
       call swap(ar%ind(1,:),ar%ind(i,:))
       call shift_down(ar,1,i-1)
    end do

  contains
    subroutine shift_down(ar,l,r)
      use my_types
      implicit none
      type(list_index)::ar
      integer,intent(in) :: l,r
      integer :: j,jold
      integer,allocatable ::a(:)

      allocate (a(size(ar%ind,2)))
      a=ar%ind(l,:)
      jold=l
      j=l+l
      do
         if (j > r) exit
         if (j < r) then
            if (ar%ind(j,1) < ar%ind(j+1,1)) then
               j=j+1
            end if
         end if
         if (a(1) >= ar%ind(j,1)) exit
         ar%ind(jold,:)=ar%ind(j,:)
         jold=j
         j=j+j
      end do
      ar%ind(jold,:)=a
    end subroutine shift_down

    subroutine swap(a,b)
      use my_types
      implicit none
      integer,dimension(:) :: a,b
      integer,allocatable ::dum (:)

      allocate(dum(size(a)))

      dum=a
      a=b
      b=dum
    end subroutine swap

  end subroutine heap_sort_list_index

  !-----------------------------------------------------------------------------!

  subroutine locate_in_list_index(xx,x,j)
    ! Given a list index xx and a value x, returns a value j such that
    ! x is between x(j) and x(j+1)
    ! based on numerical recepies in fortran
    use my_types
    implicit none

    ! subroutine arguments
    type(list_index),intent(in):: xx
    integer,intent(in) :: x
    integer,intent(out) ::j
    integer :: jl,jm,ju
    logical :: ascnd

    ascnd = (xx%ind(xx%n,1) >= xx%ind(1,1))
    jl=0
    ju=xx%n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx%ind(jm,1))) then
          jl=jm
       else
          ju=jm
       end if
    end do
    if (x == xx%ind(1,1)) then
       j=1
    else if (x == xx%ind(xx%n,1)) then
       j=xx%n
    else
       j=jl
    end if


    if (j==0) then
       j=-1
    else if (x/=xx%ind(j,1)) then
       j=-1
    end if

  end subroutine locate_in_list_index

  !-----------------------------------------------------------------------------!

  function is_in_list_index(xx,x) result(bol)
    ! Test if a given value x lies in the list xx
    use my_types
    implicit none

    ! subroutine arguments
    type(list_index) :: xx
    integer          :: x
    logical          :: bol

    ! local variables

    integer :: j

    ! maybe needed
    if (xx%n==0) then
       bol=.false.
       return
    end if

    call locate_in_list_index(xx,x,j)

    if (j==-1) then
       bol=.false.
    else
       bol=.true.
    end if

  end function is_in_list_index

  !-----------------------------------------------------------------------------!

  subroutine linear_interp_val(gr,x,y,intp)
    ! performs a 2d linear interpolation for the grid value at the position x,y
    use my_types
    implicit none
    ! subroutine arguments
    type(scalar_field_2d) :: gr
    real(kind=dbl),intent(in)  :: x,y
    real(kind=dbl),intent(out) :: intp

    !local variables
    real(kind=dbl)::bs(2),bt(2),s,t
    integer :: i,j,p,q,m,n

    p=int((x-gr%x0)/gr%dx)+1+gr%cn
    q=int((y-gr%y0)/gr%dy)+1+gr%cn

    s=(x-((p-1-gr%cn)*gr%dx+gr%x0))/gr%dx
    t=(y-((q-1-gr%cn)*gr%dy+gr%y0))/gr%dy

    bs(1)=1.0-s
    bs(2)=s

    bt(1)=1.0-t
    bt(2)=t

    intp=0.0_dbl

    do i=1,2
       do j=1,2

          m=p-1+i
          n=q-1+j

          intp=intp+bs(i)*bt(j)*gr%val(m,n)
       end do
    end do

  end subroutine linear_interp_val

  !-----------------------------------------------------------------------------!
  subroutine cubic_bspline_basis(gr,x,y,bs,bt)
    ! computes the basis functions  for the point x,y
    use my_types
    implicit none

    ! subroutine arguments
    type(scalar_field_2d) :: gr
    real(kind=dbl),intent(in)  :: x,y

    ! local variables
    real(kind=dbl)::bs(4),bt(4),s,t
    integer :: p,q

    p=int((x-gr%x0)/gr%dx)+1+gr%cn
    q=int((y-gr%y0)/gr%dy)+1+gr%cn

    s=(x-((p-1-gr%cn)*gr%dx+gr%x0))/gr%dx
    t=(y-((q-1-gr%cn)*gr%dy+gr%y0))/gr%dy

    bs(1)=1.0/6.0*(1.0-s)**3
    bs(2)=1.0/6.0*(4.0-6.0*s**2+3.0*s**3)
    bs(3)=1.0/6.0*(1.0+3.0*s+3.0*s**2-3.0*s**3)
    bs(4)=1.0/6.0*s**3

    bt(1)=1.0/6.0*(1.0-t)**3
    bt(2)=1.0/6.0*(4.0-6.0*t**2+3.0*t**3)
    bt(3)=1.0/6.0*(1+3.0*t+3.0*t**2-3.0*t**3)
    bt(4)=1.0/6.0*t**3

  end subroutine cubic_bspline_basis

  !-----------------------------------------------------------------------------!

  subroutine cubic_bspline_interp_val(gr,x,y,intp)
    !performes a 2d cubic bspline interpolation for the grid value at the position x,y
    use my_types
    implicit none

    ! subroutine arguments
    type(scalar_field_2d) :: gr
    real(kind=dbl),intent(in)  :: x,y
    real(kind=dbl),intent(out) :: intp

    ! local variables
    real(kind=dbl)::bs(4),bt(4),s,t
    integer :: p,q,m,n,i,j

    p=int((x-gr%x0)/gr%dx)+1+gr%cn
    q=int((y-gr%y0)/gr%dy)+1+gr%cn

    s=(x-((p-1-gr%cn)*gr%dx+gr%x0))/gr%dx
    t=(y-((q-1-gr%cn)*gr%dy+gr%y0))/gr%dy

    bs(1)=1.0/6.0*(1.0-s)**3
    bs(2)=1.0/6.0*(4.0-6.0*s**2+3.0*s**3)
    bs(3)=1.0/6.0*(1.0+3.0*s+3.0*s**2-3.0*s**3)
    bs(4)=1.0/6.0*s**3

    bt(1)=1.0/6.0*(1.0-t)**3
    bt(2)=1.0/6.0*(4.0-6.0*t**2+3.0*t**3)
    bt(3)=1.0/6.0*(1+3.0*t+3.0*t**2-3.0*t**3)
    bt(4)=1.0/6.0*t**3

    intp=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          intp=intp+bs(i)*bt(j)*gr%val(m,n)
       end do
    end do

  end subroutine cubic_bspline_interp_val

  !-----------------------------------------------------------------------------!

  subroutine cubic_bspline_interp_val_1der(gr,x,y,f,fx,fy)
    ! Does a cubic bspline interpolation in 2d and provides also the derviatives
    ! in the two spatial directions

    use my_types
    implicit none

    !subroutine arguments
    type(scalar_field_2d),intent(in) :: gr
    real(kind=dbl),intent(in)  :: x,y
    real(kind=dbl),intent(out) :: f,fx,fy

    ! local variables
    real(kind=dbl)::bs0(4),bt0(4),bs1(4),bt1(4),s,t
    integer :: p,q,m,n,i,j

    p=int((x-gr%x0)/gr%dx)+1+gr%cn
    q=int((y-gr%y0)/gr%dy)+1+gr%cn

    s=(x-((p-1-gr%cn)*gr%dx+gr%x0))/gr%dx
    t=(y-((q-1-gr%cn)*gr%dy+gr%y0))/gr%dy

    bs0(1)=1.0/6.0*(1.0-s)**3
    bs0(2)=1.0/6.0*(4.0-6.0*s**2+3.0*s**3)
    bs0(3)=1.0/6.0*(1.0+3.0*s+3.0*s**2-3.0*s**3)
    bs0(4)=1.0/6.0*s**3

    bt0(1)=1.0/6.0*(1.0-t)**3
    bt0(2)=1.0/6.0*(4.0-6.0*t**2+3.0*t**3)
    bt0(3)=1.0/6.0*(1+3.0*t+3.0*t**2-3.0*t**3)
    bt0(4)=1.0/6.0*t**3

    bs1(1)=-1.0/2.0*(1-s)**2
    bs1(2)=+1.0/2.0*(-4*s+3*s**2)
    bs1(3)=+1.0/2.0*(1+2*s-3*s**2)
    bs1(4)=+1.0/2.0*s**2

    bt1(1)=-1.0/2.0*(1-t)**2
    bt1(2)=+1.0/2.0*(-4*t+3*t**2)
    bt1(3)=+1.0/2.0*(1+2*t-3*t**2)
    bt1(4)=+1.0/2.0*t**2

    f=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          f=f+bs0(i)*bt0(j)*gr%val(m,n)
       end do
    end do

    fx=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          fx=fx+bs1(i)*bt0(j)*gr%val(m,n)
       end do
    end do

    fx=fx*1.0/gr%dx

    fy=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          fy=fy+bs0(i)*bt1(j)*gr%val(m,n)
       end do
    end do

    fy=fy*1.0/gr%dy
  end subroutine cubic_bspline_interp_val_1der

  !-------------------------------------------------------------------------!

  subroutine cubic_bspline_interp_val_12der(gr,x,y,f,fx,fy,fxy,fxx,fyy)
    ! Does a cubic bspline interpolation in 2d and provides also the first and
    ! second order derviatives in the two spatial directions
    use my_types
    implicit none

    !subroutine arguments
    type(scalar_field_2d),intent(in) :: gr
    real(kind=dbl),intent(in)  :: x,y
    real(kind=dbl),intent(out) :: f,fx,fy,fxx,fyy,fxy

    ! local variables
    real(kind=dbl)::bs0(4),bt0(4),bs1(4),bt1(4),bs2(4),bt2(4),s,t
    integer :: p,q,m,n,i,j


    p=int((x-gr%x0)/gr%dx)+1+gr%cn
    q=int((y-gr%y0)/gr%dy)+1+gr%cn

    s=(x-((p-1-gr%cn)*gr%dx+gr%x0))/gr%dx
    t=(y-((q-1-gr%cn)*gr%dy+gr%y0))/gr%dy

    bs0(1)=1.0/6.0*(1.0-s)**3
    bs0(2)=1.0/6.0*(4.0-6.0*s**2+3.0*s**3)
    bs0(3)=1.0/6.0*(1.0+3.0*s+3.0*s**2-3.0*s**3)
    bs0(4)=1.0/6.0*s**3

    bt0(1)=1.0/6.0*(1.0-t)**3
    bt0(2)=1.0/6.0*(4.0-6.0*t**2+3.0*t**3)
    bt0(3)=1.0/6.0*(1+3.0*t+3.0*t**2-3.0*t**3)
    bt0(4)=1.0/6.0*t**3

    bs1(1)=-1.0/2.0*(1.0-s)**2
    bs1(2)=+1.0/2.0*(-4.0*s+3.0*s**2)
    bs1(3)=+1.0/2.0*(1.0+2.0*s-3.0*s**2)
    bs1(4)=+1.0/2.0*s**2

    bt1(1)=-1.0/2.0*(1.0-t)**2
    bt1(2)=+1.0/2.0*(-4.0*t+3.0*t**2)
    bt1(3)=+1.0/2.0*(1.0+2.0*t-3.0*t**2)
    bt1(4)=+1.0/2.0*t**2

    bs2(1)=+(1-s)
    bs2(2)=+1.0/2.0*(-4.0+6.0*s)
    bs2(3)=+1.0/2.0*(2.0-6.0*s)
    bs2(4)=+1.0*s

    bt2(1)=+(1-t)
    bt2(2)=+1.0/2.0*(-4.0+6.0*t)
    bt2(3)=+1.0/2.0*(2.0-6.0*t)
    bt2(4)=+1.0*t

    f=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          f=f+bs0(i)*bt0(j)*gr%val(m,n)
       end do
    end do

    fx=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          fx=fx+bs1(i)*bt0(j)*gr%val(m,n)
       end do
    end do

    fx=fx*1.0/gr%dx

    fy=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          fy=fy+bs0(i)*bt1(j)*gr%val(m,n)
       end do
    end do

    fy=fy*1.0/gr%dy

    fxx=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          fxx=fxx+bs2(i)*bt0(j)*gr%val(m,n)
       end do
    end do

    fxx=fxx*1.0/gr%dx**2

    fyy=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          fyy=fyy+bs0(i)*bt2(j)*gr%val(m,n)
       end do
    end do

    fyy=fyy*1.0/gr%dy**2


    fxy=0.0_dbl

    do i=1,4
       do j=1,4
          m=p-2+i
          n=q-2+j
          if (m>gr%nx+2*gr%cn) then
             m=gr%nx
          end if
          if (n>gr%ny+2*gr%cn) then
             n=gr%ny
          end if
          if (m<1) then
             m=1
          end if
          if (n<1) then
             n=1
          end if
          fxy=fxy+bs1(i)*bt1(j)*gr%val(m,n)
       end do
    end do

    fxy=fxy*1.0/(gr%dx*gr%dy)

  end subroutine cubic_bspline_interp_val_12der

  !------------------------------------------------------------------------------!

subroutine upsample_bspline(ifd,ofd,nrx,nry)

use my_types
implicit none

type(scalar_field_2d):: ifd
type(scalar_field_2d):: ofd
integer :: nrx,nry

real(kind=dbl):: x,y
integer :: i,j

  ofd%x0=ifd%x0;ofd%y0=ifd%y0
  ofd%nx=nrx*(ifd%nx-1)+1;ofd%ny=nry*(ifd%ny-1)+1
  ofd%dx=((ifd%nx-1)*ifd%dx)/(ofd%nx-1);ofd%dy=((ifd%ny-1)*ifd%dy)/(ofd%ny-1)
  ofd%cn=ifd%cn
  allocate(ofd%val(ofd%nx+ofd%cn*2,ofd%ny+ofd%cn*2))


  do i=1+ofd%cn,ofd%nx+ofd%cn
     do j=1+ofd%cn,ofd%ny+ofd%cn
        x=(i-1-ofd%cn)*ofd%dx+ofd%x0
        y=(j-1-ofd%cn)*ofd%dy+ofd%y0
        call cubic_bspline_interp_val(ifd,x,y,ofd%val(i,j))
     end do
  end do

end subroutine upsample_bspline
!--------------------------------------------------------------------------------!

  subroutine velnod2id(i,j,vmod,id)
    use my_types
    implicit none
    type(scalar_field_2d):: vmod
    integer :: i,j,id

    id=(i-1)*(vmod%ny+vmod%cn*2)+j

  end subroutine velnod2id

  !----------------------------------------------------------------------------!

  subroutine id2velnod(id,vmod,i,j)
    use my_types
    implicit none
    type(scalar_field_2d):: vmod
    integer :: i,j,id

    i=int((id-1)/(vmod%ny+vmod%cn*2))+1
    j=id-(i-1)*(vmod%ny+vmod%cn*2)

  end subroutine id2velnod

  !----------------------------------------------------------------------------!

  subroutine id2nei(id,vmod,nei)

    use my_types
    implicit none

    ! subroutine arguments
    integer :: id
    type(scalar_field_2d):: vmod
    type(list_index):: nei

    ! local variables
    integer :: near(2,4)
    integer :: i,j
    integer :: h,k
    integer :: ii,jj

    ! initialise the neighbour index
    near(1,:)=(/-1,0,+1,0/)
    near(2,:)=(/0,-1,0,+1/)

    nei%n=0

    ! localise the node in the model
    call id2velnod(id,vmod,i,j)

    do h=1,4 ! loop over the four neighbours
       ii=i+near(1,h)
       jj=j+near(2,h)
       if (ii>=1.and.ii<=vmod%nx+vmod%cn*2.and.&
            jj>=1.and.jj<=vmod%ny+vmod%cn*2) then
          nei%n=nei%n+1
          call velnod2id(ii,jj,vmod,nei%ind(nei%n,1))
       end if
    end do

  end subroutine id2nei

  !----------------------------------------------------------------------------!

  subroutine raypath_linear_interp_3d(t1,x1,y1,z1,t2,x2,y2,z2,tt,xx,yy,zz)
    ! interpolates a 3d path for a given time tt
    use my_types
    implicit none
    ! subroutine arguments
    real(kind=dbl) t1,x1,y1,z1,t2,x2,y2,z2,tt,xx,yy,zz

    xx=linear_interp_1d(t1,x1,t2,x2,tt)
    yy=linear_interp_1d(t1,y1,t2,y2,tt)
    zz=linear_interp_1d(t1,z1,t2,z2,tt)

  end subroutine raypath_linear_interp_3d

  !--------------------------------------------------------------------------!

  function linear_interp_1d (x1,y1,x2,y2,x) result(y)
    ! computes the corresponding y value for  the input parmeter x using a
    ! linear interpolation based on the points (x1,y1) and (x2,y2)
    use my_types
    implicit none
    !function arguments
    real(kind=dbl),intent(in):: x1,x2,y1,y2,x
    real(kind=dbl)::y
    !local variables
    real(kind=dbl)::a,b
    a=(y2-y1)/(x2-x1)
    b=(y1-a*x1)
    y=a*x+b
  end function linear_interp_1d

  !--------------------------------------------------------------------------!

  function in_polygon(x,y,x0,y0) result(ans)
    ! this function test if a given node lies in a 2d polygon. ans=1 if the point
    ! lies on the side or in the polygon otherwise ans=0. This subroutine
    ! calls a much more complex subroutine
    use my_types
    implicit none
    ! subroutine arguments
    real(kind=dbl)::x(:),y(:),x0,y0
    integer :: ans

    !local variables
    integer :: m,n
    real(kind=dbl),allocatable :: xx(:),yy(:)



    if (x0>maxval(x).or.x0<minval(x).or.y0>maxval(y).or.y0<minval(y)) then
       ans=-1
       return
    end if


    n=size(x)
    allocate(xx(n+1))
    allocate(yy(n+1))
    xx(1:n)=x;xx(n+1)=x(1)
    yy(1:n)=y;yy(n+1)=y(1)

    call point_polygon_relation (xx,yy,x0,y0,ans,m)

  end function in_polygon

  !---------------------------------------------------------------------------!

  subroutine point_polygon_relation (x, y, x0, y0, k, m)
    ! given a polygonal line connecting the vertices (x(:)y(:)) tkaen in the order
    ! from 1 to n. it is assumes that he polygonal path is a loop where x(n)=x(1)
    ! and y(n)=y(1) or there is an arc from (x(n),y(n)) to (x(1),y(1)).
    ! the polygon may cross itself any number of times
    ! (x0,y0) is an arbitrary point l and m are variables
    ! On output, l and m are assigned the following values...
    !    k = -1   if (x0,y0) is out side the polygonal path
    !    k =  0   if (x0,y0) lies on the polygonal path
    !    k =  1   if (x0,y0) lie in side the polygonal path
    ! m=0 if (x0,y0) is on or outside the paths. if (x0,y0) is inside the
    ! path then m is the wnding number of the path around the point.

    use my_types
    use my_constants

    implicit none
    ! subroutine arguments
    real(kind=dbl), intent(in) :: x0, y0, x(:), y(:)
    integer,intent(out)        :: k, m

    ! local variabels
    integer:: i, n
    real(kind=dbl)    :: angle, eps, pi2, sum, theta, theta1, thetai, tol, u, v

    ! eps is a machin dependent constant, eps is the smallest number so
    ! that 1.0 + EPS > 1.0

    eps = EPSILON(1.0_dbl)

    pi2 = 2.0*pi           ! 2 pi just for making the life easier
    tol = 4.0*eps*pi       ! angle tolerance
    k = -1                 ! initialize as outside
    m = 0                  ! outside of the paths

    n=size(x)              ! number of corners of the polygon

    if (x(1)==x(n).and.y(1)==y(n))then
       n=n-1
    end if

    ! check if x0,y0 coincides with first corner of the polygon
    u = x(1)-x0;v = y(1)-y0
    If (u == 0.0 .and. v==0.0) then
       k=0;return
    end if

    theta1 = datan2(v, u)

    sum = 0.0
    theta = theta1
    do i = 2, n ! loop over all the vertices
       ! distance from point to current polygon point
       u = x(i) - x0
       v = y(i) - y0
       ! check if x0,y0 coincides with the current corner of the polygon
       if (u==0.0 .and. v==0.0) then
          k=0 ; return
       end if
       ! update angle
       thetai = datan2(v, u)
       angle = abs(thetai - theta)
       ! check if point lies on vertice
       if (abs(angle - pi) < tol) then
          return
       end if
       ! if angle > pi reset to corresponding smaller one
       if (angle > pi)  then
          angle = angle - pi2
       end if
       ! check if we are going counter clockwise
       if (theta > thetai) then
          angle = -angle
       end if
       ! update the angle sum which is used to count the windings
       sum = sum + angle
       theta = thetai
    end do

    !check between last and first node of the polygon
    angle = abs(theta1 - theta)
    if (abs(angle - pi) < tol) then
       return
    end if
    ! if angle > pi reset to corresponding smaller one
    if (angle > pi) then
       angle = angle - pi2
    end if
    ! check if we are going counter clockwise
    if (theta > theta1) then
       angle = -angle
    end if
    ! update the angle sum which is used to count the windings
    sum = sum + angle
    ! sum =2*pi*m where m is the winding number
    m = abs(sum)/pi2 + 0.2
    if (m == 0) then
       return
    end if
    k = 1  !the point is inside the polygon
    ! check if we were walking counter clockwise
    ! if yes change the sign of the winding number
    if (sum < 0.0_dbl) then
       m = -m;return
    end if

  end subroutine point_polygon_relation

  !---------------------------------------------------------------------------!

  subroutine closest_point_on_line(x,y,x1,y1,x2,y2,xcp,ycp)
    ! Given a line between x1,y1 and x2,y2 this subroutien calculates the
    ! closest point on that line to the point x,y.
    use my_types
    implicit none

    ! subroutine arguments
    real(kind=dbl) :: x,y,x1,y1,x2,y2
    real(kind=dbl) :: xcp,ycp

    ! local variables
    real(kind=dbl) :: a1,b1,a2,b2

    ! first check if line is vertical or horizontal
    if (y1==y2) then
       xcp=x
       ycp=y1
    else if(x1==x2) then
       xcp=x1
       ycp=y
    else
       a1=(y2-y1)/(x2-x1)
       b1=y1-(x1*a1)
       a2=-1.0_dbl/a1
       b2=y-(x*a2)
       xcp=-(b1-b2)/(a1-a2)
       ycp=-(a2*b1-a1*b2)/(a1-a2)
    end if

  end subroutine closest_point_on_line

  !---------------------------------------------------------------------------!

  function great_circle_distance(x1,x2,y1,y2) result(d)

    use my_types
    use my_constants
    implicit none
    real(kind=dbl):: x1,x2,y1,y2
    real(kind=dbl):: x1r,x2r,y1r,y2r
    real(kind=dbl):: dx,dy
    real(kind=dbl):: ds,d


    x1r=x1/radian
    x2r=x2/radian

    y1r=y1/radian
    y2r=y2/radian

!    ds=acos(cos(y1r)*cos(y2r)*cos(x1r-x2r)+sin(y1r)*sin(y2r))
!write(*,*) ds
    ds=2.0*asin(sqrt(sin((y2r-y1r)/2.0)**2+cos(y1r)*cos(y2r)*sin((x2r-x1r)/2.0)**2))
!write(*,*) ds
!pause
    d=rearth*ds

  end function great_circle_distance

  !---------------------------------------------------------------------------!
  subroutine normalize_vector(v)
    use my_types
    implicit none
    ! subroutine arguments
    real(kind=dbl):: v(:)
    ! local variables
    real(kind=dbl):: s
    integer :: i

    s=0.0_dbl

    do i=1,ubound(v,1)
       s=s+v(i)**2
    end do

    s=sqrt(s)

    do i=1,ubound(v,1)
       v(i)=v(i)/s
    end do

  end subroutine normalize_vector

  function cross_product(a,b) result(c)
    ! cross product
    use my_types
    implicit none
    ! subroutine arguments
    real(kind=dbl):: a(3),b(3)
    real(kind=dbl) :: c(3)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(1)*b(3)-a(3)*b(1)
    c(3)=a(1)*b(2)-a(2)*b(1)
  end function cross_product

  !---------------------------------------------------------------------------!

  function dot_product(a,b) result(dtp)
    ! dot product for two n dimensional vectors
    use my_types
    implicit none
    ! subrotuine arguments
    real(kind=dbl),intent(in)::a(:),b(:)
    real(kind=dbl) dtp
    ! local variables
    integer :: i
    dtp=0.0_dbl
    do i=1,size(a)
       dtp=dtp+a(i)*b(i)
    end do
  end function dot_product

  !--------------------------------------------------------------------------!

  function rms_error(d1,d2) result(rms)
    ! computes the rms error between d1 and d2
    use my_types
    implicit none
    !subroutine arguments
    real(kind=dbl) :: d1(:),d2(:)
    real(kind=dbl) :: rms
    ! local variables
    integer :: i,n

    n=size(d1,1)
    rms=0.0_dbl
    do i=1,n
       rms=rms+(d1(i)-d2(i))**2
    end do

    rms=sqrt(rms/n)

  end function rms_error

  !--------------------------------------------------------------------------!

  real(kind=dbl) function gasdev(idummy)
    ! Returns a normally distributed deviate with zero mean and unit variance,
    ! using ran1(idum)  as thesource of uniformdeviates.
    use my_types
    implicit none
    integer :: iset,i,idummy
    integer, parameter :: imax=100000
    real(kind=dbl) :: fac,rsq,v1,v2
    real(kind=dbl), SAVE :: gset
    iset=0
    if(iset==0)then
       do i=1,imax
          v1=2.0*ran1(idummy)-1.
          v2=2.0*ran1(idummy)-1.
          rsq=v1**2+v2**2
          if(rsq<1.and.rsq/=0.) exit
       end do
       fac=sqrt(-2.0*log(rsq)/rsq)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
    else
       gasdev=gset
       iset=0
    end if
  end function gasdev

  !--------------------------------------------------------------------------!

  real(kind=dbl) function ran1(iseed)
    ! Minimal randomnumber generator of Parka nd Miller with Bays-Durham shuffle and
    ! addedsafeguards. Returns auniformrandomdeviate between0.0and1.0(exclusive of
    ! the endpoint values). Call with idum a negative integer to initialize; thereafter,
    ! donot  alter idum between successive deviates in a sequence. RNMX should
    ! approximate the largest floating value that is less than 1.
    use my_types
    implicit none
    integer :: iseed,j,k
    integer, parameter :: ia=16807,im=2147483647,iq=127773
    integer, parameter :: ir=2836,ntab=32,ndiv=1+(im-1)/ntab
    integer, save :: iy,iv(ntab)
    real(kind=dbl), parameter :: eps=1.2e-7,rnmx=1.0-eps,am=1./im
    iv=ntab*0
    iy=0
    if(iseed<=0.or.iy==0)then
       do j=ntab+8,1,-1
          k=iseed/iq
          iseed=ia*(iseed-k*iq)-ir*k
          if(iseed<0)iseed=iseed+im
          if(j<=ntab)iv(j)=iseed
       end do
       iy=iv(1)
    end if
    k=iseed/iq
    iseed=ia*(iseed-k*iq)-ir*k
    if(iseed<0)iseed=iseed+im
    j=1+iy/ndiv
    iy=iv(j)
    iv(j)=iseed
    ran1=min(am*iy,rnmx)
  end function ran1

!----------------------------------------------------------------------------

  real(kind=dbl)  FUNCTION ran2(iseed)

    use my_types
    implicit none

    INTEGER iseed,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL(kind=dbl):: AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
         NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER iseed2,j,k,iv(NTAB),iy
    SAVE iv,iy,iseed2
    DATA iseed2/123456789/, iv/NTAB*0/, iy/0/
    if (iseed.le.0) then
       iseed=max(-iseed,1)
       iseed2=iseed
       do  j=NTAB+8,1,-1
          k=iseed/IQ1
          iseed=IA1*(iseed-k*IQ1)-k*IR1
          if (iseed.lt.0) iseed=iseed+IM1
          if (j.le.NTAB) iv(j)=iseed
       end do
       iy=iv(1)
    endif
    k=iseed/IQ1
    iseed=IA1*(iseed-k*IQ1)-k*IR1
    if (iseed.lt.0) iseed=iseed+IM1
    k=iseed2/IQ2
    iseed2=IA2*(iseed2-k*IQ2)-k*IR2
    if (iseed2.lt.0) iseed2=iseed2+IM2
    j=1+iy/NDIV
    iy=iv(j)-iseed2
    iv(j)=iseed
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
  end FUNCTION ran2


  !------------------------------------------------------------------------!

  subroutine pearsn(x,y,r,prob,z)
    ! Given two arrays x and y thsi subroutine computes their correlation
    ! coeffcien R, the significance leve at which the null hypothesis of zero
    ! correlation is disproved (prob whos small values indicate a signifcant
    ! correlation), and Fisher's z, whos calue can be used in furter statistical
    ! tests
    use my_types

    implicit none

    ! subroutine arguments
    real(kind=dbl) :: x(:),y(:)
    real(kind=dbl) :: r,prob,z

    ! local variables
    real(kind=dbl),parameter:: tiny=1.e-20
    integer n,j
    real(kind=dbl):: ax,ay,df,sxx,sxy,syy,t,xt,yt
    n=size(x)


    ! comput the mean
    ax=0.01_dbl; ay=0.0_dbl
    do  j=1,n
       ax=ax+x(j)
       ay=ay+y(j)
    end do
    ax=ax/n
    ay=ay/n

    ! compute he correlation coefficient
    sxx=0.0_dbl; syy=0.0_dbl;sxy=0.0_dbl
    do  j=1,n
       xt=x(j)-ax
       yt=y(j)-ay
       sxx=sxx+xt**2
       syy=syy+yt**2
       sxy=sxy+xt*yt
    end do

    r=sxy/(sqrt(sxx*syy)+tiny)
    z=0.5*log(((1.+r)+tiny)/((1.-r)+tiny))
    df=n-2
    t=r*sqrt(df/(((1.-r)+tiny)*((1.+r)+tiny)))
    prob=betai(0.5*df,0.5_dbl,df/(df+t**2))
    return
  end subroutine pearsn

  !--------------------------------------------------------------------------!

  real(kind=dbl) function betai(a,b,x)
    !returns the incomplete beta functions /x(a,B)
    use my_types
    implicit none
    ! fucntion arguments
    real(kind=dbl):: a,b,x
    ! local variables
    real(kind=dbl)::  bt
    if(x<0.0_dbl.or.x>1.0_dbl) then
       write(*,*) 'bad argument x in betai'
    end if
    if(x==0.0_dbl.or.x==1.0_dbl)then
       bt=0.
    else ! factors in fornt of the continued fraction
       bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
    endif
    if(x.lt.(a+1.)/(a+b+2.))then  ! use continued fraction directly
       betai=bt*betacf(a,b,x)/a
       return
    else
       betai=1.-bt*betacf(b,a,1.-x)/b ! use continued fraction after making the symmetry transformation
       return
    endif
  end function betai

  !--------------------------------------------------------------------------!

  real(kind=dbl)  function betacf(a,b,x)
    !    continued fraction for incomplete beta function, used by betai
    use my_types
    implicit none
    ! function arguments
    real(kind=dbl)::a,b,x
    ! local variables
    integer,parameter :: maxit=100
    real(kind=dbl),parameter :: eps=3.e-7,fpmin=1.e-30
    integer m,m2
    real(kind=dbl) :: aa,c,d,del,h,qab,qam,qap

    !These Q's will be use in factors which occur in the coefficients
    qab=a+b
    qap=a+1.
    qam=a-1.
    ! First step of Lentz's method
    c=1.
    d=1.-qab*x/qap
    if(abs(d).lt.fpmin)d=fpmin
    d=1./d
    h=d
    do m=1,maxit
       ! continued fraction evaluation by the recurrence method
       m2=2*m
       aa=m*(b-m)*x/((qam+m2)*(a+m2))
       d=1.+aa*d
       ! one step (the even one) of the recurrence
       if(abs(d).lt.fpmin)d=fpmin
       c=1.+aa/c
       if(abs(c).lt.fpmin)c=fpmin
       d=1./d
       h=h*d*c
       aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
       d=1.+aa*d
       ! next steop of the recurrence (the odd one)
       if(abs(d).lt.fpmin)d=fpmin
       c=1.+aa/c
       if(abs(c).lt.fpmin)c=fpmin
       d=1./d
       del=d*c
       h=h*del
       if(abs(del-1.).lt.eps) then
          ! ar we done ?
          betacf=h
          return
       end if
    end do
    write(*,*) 'a or b too big, or maxit too small in betacf'
    betacf=h

  end function betacf

  !--------------------------------------------------------------------------!

  real(kind=dbl) function gammln(xx)
    ! Returns the value ln(Gammma(xx)) for xx>0
    use my_types
    implicit none

    real(kind=dbl)::xx
    integer:: j
    real(kind=dbl):: ser,tmp,x,y

    real(kind=dbl) :: stp = 2.5066282746310005_dbl
    real(kind=dbl) :: cof(6) = (/76.18009172947146_dbl,&
         -86.50532032941677_dbl,24.01409824083091_dbl,&
         -1.231739572450155_dbl,0.1208650973866179e-2_dbl,&
         -0.5395239384953e-5_dbl/)
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do  j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
    end do
    gammln=tmp+log(stp*ser/x)
  end function  gammln








!-----------------------------------------------------------------------------!


   subroutine moment(data,n,ave,adev,sdev,var,skew,curt)
      INTEGER n
      real(kind=dbl):: adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      real(kind=dbl):: p,s,ep
      s=0.
      do j=1,n
        s=s+data(j)
end do
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
end do
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      endif
      return
      end subroutine

!------------------------------------------------------------------------------!














end module my_functions

