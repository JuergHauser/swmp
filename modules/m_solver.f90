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
!   along with swmp. If not, see <http://www.gnu.org/licenses/>.
!
!   You may contact the author by:
!       e-mail:  juerg.hauser@csiro.au
!
! * m_solver *
!
! solvers for the kineamtic ray tracing equations on a spherical surface
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


module m_solver

contains


  subroutine derivs_k(x,dxdt,vmod)
    ! odes for the ray tracing in phase space
    use my_types
    use my_functions
    use my_constants
    implicit none
    !subroutine arguments
    real(kind=dbl)::x(3),dxdt(3)
    type(scalar_field_2d) :: vmod
    !local variables
    real(kind=dbl) :: c,cx,cy
    real(kind=dbl):: colat

    call cubic_bspline_interp_val_1der(vmod,x(1),x(2),c,cx,cy)

    ! dxdt(1)=c*dcos(x(3))
    ! dxdt(2)=c*dsin(x(3))
    ! dxdt(3)=cx*dsin(x(3))-cy*dcos(x(3))
    
    ! colatitude
    ! write(*,*) c,cx,cy

cx=cx*radian
cy=cy*radian

! we are using colatitude
if (x(2)>0.0) then
   cy=-cy
end if


    ! cx=0.0_dbl
    ! cy=0.0_dbl
    ! colat=sign((90-abs(x(2)))/radian,x(2))
    colat=(90.0_dbl-abs(x(2)))/radian
 
    !colat=pi/2.0-x(2)/radian      
        
  !  dxdt(1)=radian*c/(rearth*sin(colat))*sin(x(3))
  !  dxdt(2)=radian*c/rearth*cos(x(3))
  !  dxdt(3)=sin(x(3))/rearth*cy-cos(x(3))/(rearth*sin(colat))*cx - c/rearth*sin(x(3))*1.0_dbl/tan(colat)
    

dxdt(1)=c*sin(x(3))/(rearth*sin(colat))*radian
dxdt(2)=-c*cos(x(3))/rearth*radian
dxdt(3)=-1.0/rearth*(sin(x(3))*cy+cos(x(3))/sin(colat)*cx)+c/rearth*sin(x(3))/tan(colat)



! we are using colatitude
if (x(2)>0.0) then
   dxdt(2)=-dxdt(2)
end if






  !  write(*,*) dxdt(3)
    
  end subroutine derivs_k
  !--------------------------------------------------------------------------------!

  subroutine rk1_k(x,xout,h,vmod)
    ! first order runge kutta time step. This is the same as a single eulerian time
    ! step

    use my_types
    implicit none

    ! subroutine arguments
    real(kind=dbl),intent(in) :: x(3),h
    type(scalar_field_2d),intent(in) :: vmod
    real(kind=dbl),intent(out)::xout(3)
    ! local variables
    real(kind=dbl) :: dxdt(3)


    call derivs_k(x,dxdt,vmod)
    xout=x+h*dxdt

  end subroutine rk1_k
  !---------------------------------------------------------------------------------!
  subroutine rk4_k(x,xout,h,vmod)	
    ! given values for n variables Y and their derivatives dydx known at x, use the forth
    ! order Rugne Kutta method to advance the solution of an interval h and return the
    ! incremented variables as xout. The user supplies the subroutine derivs(x,ydxdy) 
    ! which returns deivatives dydx at X.	

    use my_types
    implicit none

    !subroutine arguments
    real(kind=dbl),intent(in):: x(3),h
    type(scalar_field_2d),intent(in) :: vmod
    real(kind=dbl),intent(out)::xout(3)

    !local variables
    real(kind=dbl):: h6,hh,dxdt(3),dxt(3),dxm(3),xt(3)

    hh=h*0.5_dbl
    h6=h/6.0_dbl

    call derivs_k(x,dxdt,vmod) ! first step
    xt=x+hh*dxdt
    call derivs_k(xt,dxt,vmod)  ! second step
    xt=x+hh*dxt              
    call derivs_k(xt,dxm,vmod)  ! third step
    xt=x+h*dxm
    dxm=dxt+dxm
    call derivs_k(xt,dxt,vmod)  ! fourth step
    xout=x+h6*(dxdt+dxt+2.0_dbl*dxm)
  end subroutine rk4_k

  !------------------------------------------------------------------------------------!

  subroutine rkck5_k(x,xout,xerr,h,vmod)
    ! Given values for x(n) and the subroutine derivs which provides dxdt(n), us the!
    ! fifth order Cash-Karp Runge Kutta Method do advance the solution over an interval
    ! h and returns the incremented varaibles a xout. An estimation of the locla truncation
    ! error in xout is also provided using the embedded fourth order method. 
    use my_types
    implicit none
    ! subroutien arguments
    type(scalar_field_2d),intent(in)::vmod
    real(kind=dbl),intent(in) ::x(3),h
    real(kind=dbl),intent(out)::xout(3),xerr(3)

    !local variables
    real(kind=dbl) :: dxdt(3),ak2(3),ak3(3),ak4(3),ak5(3),ak6(3),xtemp(3)
    real(kind=dbl), parameter :: A2=0.2_dbl,A3=0.3_dbl,A4=0.6_dbl,A5=1.0_dbl,&
         A6=0.875_dbl,B21=0.2_dbl,B31=3.0_dbl/40.0_dbl,B32=9.0_dbl/40.0_dbl,&
         B41=0.3_dbl,B42=-0.9_dbl,B43=1.2_dbl,B51=-11.0_dbl/54.0_dbl,&
         B52=2.5_dbl,B53=-70.0_dbl/27.0_dbl,B54=35.0_dbl/27.0_dbl,&
         B61=1631.0_dbl/55296.0_dbl,B62=175.0_dbl/512.0_dbl,&
         B63=575.0_dbl/13824.0_dbl,B64=44275.0_dbl/110592.0_dbl,&
         B65=253.0_dbl/4096.0_dbl,C1=37.0_dbl/378.0_dbl,&
         C3=250.0_dbl/621.0_dbl,C4=125.0_dbl/594.0_dbl,&
         C6=512.0_dbl/1771.0_dbl,DC1=C1-2825.0_dbl/27648.0_dbl,&
         DC3=C3-18575.0_dbl/48384.0_dbl,DC4=C4-13525.0_dbl/55296.0_dbl,&
         DC5=-277.0_dbl/14336.0_dbl,DC6=C6-0.25_dbl

    call derivs_k(x,dxdt,vmod)
    xtemp=x+B21*h*dxdt
    call derivs_k(xtemp,ak2,vmod)
    xtemp=x+h*(B31*dxdt+B32*ak2)
    call derivs_k(xtemp,ak3,vmod)
    xtemp=x+h*(B41*dxdt+B42*ak2+B43*ak3)
    call derivs_k(xtemp,ak4,vmod)
    xtemp=x+h*(B51*dxdt+B52*ak2+B53*ak3+B54*ak4)
    call derivs_k(xtemp,ak5,vmod)
    xtemp=x+h*(B61*dxdt+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
    call derivs_k(xtemp,ak6,vmod)
    xout=x+h*(C1*dxdt+C3*ak3+C4*ak4+C6*ak6)
    xerr=h*(DC1*dxdt+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
  end subroutine rkck5_k

  !------------------------------------------------------------------------------------!

  subroutine rkqs5_k(x,xout,eps,xscal,htry,hdid,hnext,vmod)
    ! fifth order Runge Kutta step with monitoring of local truncation error to ensure
    ! accuracy and adjust step size. Input are the dpendent varaibel vaector x(n) and the
    ! the derivatives dxdt(n) provided by the user derived function derivs. The to be 
    ! attempted stepsize is htry, the required accuracy eps, and the vector yscal(n)
    ! against  which the error is scaled. xout is the output and hdid is the stepsize
    ! that was actually accomplished, and hnext is the estimated next stepsize.

    use my_types

    implicit none
    ! subroutine arguments    
    type(scalar_field_2d),intent(in):: vmod   
    real(kind=dbl),intent(in) ::x(3),xscal(3),htry,eps
    real(kind=dbl),intent(out):: xout(3),hdid,hnext

    ! local variables
    real(kind=dbl),parameter::safety=0.9_dbl,pgrow=-0.2_dbl,pshrink=-0.25_dbl,errcon=1.89e-4_dbl 
    real(kind=dbl)::htemp,h,xtemp(3),xerr(3),errmax

    h=htry
    do
       call rkck5_k(x,xtemp,xerr,h,vmod)
       errmax=maxval(abs(xerr/xscal))/eps
       if (errmax <= 1.0) exit
       htemp=safety*h*(errmax**pshrink)
       h=sign(max(abs(htemp),0.1_dbl*abs(h)),h)
       if (h == 0) write(*,*) ('stepsize underflow in rkqs')
    end do
    if (errmax > errcon) then
       hnext=safety*h*(errmax**pgrow)
    else
       hnext=5.0_dbl*h
    end if
    hdid=h
    xout=xtemp
  end subroutine rkqs5_k

  !------------------------------------------------------------------------------------!

  subroutine odeint_k(xin,xout,dt,eps,hmin,vmod)
    ! Runge kutta drive with adaptive stepsize control. Integrate the n starting values xstart 
    ! over distance dt with accuracy eps. h1 is the guessed first step size and hmin is the
    ! minimum allowed stepsize(can be zero). The output is provided in the varaible xout. 
    use my_types
    implicit none


    ! subroutine arguments
    type(scalar_field_2d),intent(in):: vmod
    real(kind=dbl) :: xin(3),dt,eps,hmin
    real(kind=dbl),intent(out)::xout(3)

    !local variables
    real(kind=dbl):: x(3),dxdt(3),xscal(3),t,h,hdid,hnext,t1,t2
    integer:: nstp
    integer,parameter:: maxstp=10,tiny=1.0e-30_dbl

    x=xin
    h=dt
    t1=0_dbl
    t2=dt
    do nstp=1,maxstp
       call derivs_k(x,dxdt,vmod)
       xscal=abs(x)+abs(h*dxdt)+tiny
       if ((t+h-t2)*(t+h-t1)>0.0) h=t2-t        ! if step size can overshoot decrease
       call rkqs5_k(x,xout,eps,xscal,h,hdid,hnext,vmod) 
       t=t+hdid;x=xout ! updating time and position
       if ((t-t2)*(t2-t1) >= 0.0) then
          xout=x
          return !normal exit
       end if
       if (abs(hnext) < hmin) then
          write(*,*)('stepsize smaller than minimum in odeint'),hnext
       end if
       h=hnext
    end do
    write(*,*) 'too many steps in odeint'

  end subroutine odeint_k





end module m_solver



