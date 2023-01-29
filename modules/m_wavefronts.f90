! * m_wavefronts *
!
! This module contans subroutine using the Lagrangian approach for the update
! step of the bicharacteristic strip.
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

module m_wavefronts


contains

  !---------------------------------------------------------------------------------!

  subroutine delete_strip(head,tail)
    ! this subroutine erases a bicharacteristic strip if it is present
    use my_types
    implicit none
    ! subroutine arguments
    type(node),pointer :: head,tail
    !local variables
    type(node), pointer :: oldnode=> null()
    type(node),pointer  :: cur => null()

    if (associated(head)) then
       do while(head%next%pid/=tail%pid)
          cur=>head%next
          ! remove in between
          oldnode => cur
          cur%prev%next => cur%next
          cur%next%prev => cur%prev
          cur%next%con(1)=0
          cur%prev%con(2)=0
          cur =>head
          deallocate(oldnode)
       end do

       ! we now only have the head and the tail
       ! we just need to get rid of these two
       deallocate(head)
       deallocate(tail)
       nullify(head)
       nullify(tail)
    end if

  end subroutine delete_strip

  !---------------------------------------------------------------------------------!

  subroutine setup_strip_for_point_source(vmod,sous,recs,params,head,tail)
    ! this subroutine creates the bicharactersitic strip for a given source point
    ! params%sid out of the source point list given in source%pos 
    use my_types
    use my_constants
    use my_functions
    implicit none
    ! subroutine arguments
    type(rat_conf)    :: params
    type(node),pointer :: head,tail
    type(scalar_field_2d) :: vmod
   ! type(wavefronts)   :: wafs
    type(sources)      :: sous
    type(receivers)    :: recs
    !local variables
    integer             :: h,i
    real(kind=dbl)      :: x,y,z

    real(kind=dbl):: v0

    call delete_strip(head,tail)


       x=sous%pos(params%sid,1)
       y=sous%pos(params%sid,2)

call cubic_bspline_interp_val(vmod,x,y,v0)

    !rat_setup bicharacteristic strip storage structure
    h=0
    do i=1,params%wt%n_senod
       z=  -pi+(2*pi)/(params%wt%n_senod+1)*(i) 
       if (.not. associated(head)) then
          allocate(head)
          tail=>head
          nullify(tail%next)
          nullify(tail%prev)
          tail%pos(1)=x
          tail%pos(2)=y
          tail%pos(3)=z
          tail%con(1)=1
          tail%con(2)=1
          h=h+1
          tail%pid=h
          tail%spf=0
tail%dinr=0
        !  wafs%stat(h,1)=1
       else
          allocate(tail%next)
          tail%next%prev=>tail
          tail => tail%next
          nullify(tail%next)
          tail%pos(1)=x
          tail%pos(2)=y
          tail%pos(3)=z
          tail%con(1)=1
          tail%con(2)=1
          h=h+1
          tail%pid=h
        !     wafs%stat(h,1)=1
tail%spf=0
tail%dinr=0
       end if
    end do

    params%wt%n_nod=params%wt%n_senod
    params%wt%nfrid=h+1

    params%wt%pdmi=(2.0*pi)/params%wt%n_nod*0.5
    params%wt%pdma=(2.0*pi)/params%wt%n_nod*2.0

    params%wt%scax=(2.0*pi)/(vmod%nx*vmod%dx-vmod%x0)
    params%wt%scay=(2.0*pi)/(vmod%ny*vmod%dy-vmod%y0)

  end subroutine setup_strip_for_point_source

!----------------------------------------!

  subroutine reset_wavefronts_and_arrivals(wafc,recs,recmode,conf)
    ! this subroutine erases wavefronts and the receiver recordings in the variables
    ! recs and was. This variables are than reset to their initial values so that the
    ! computation for the next source can be started.
    use my_types
    implicit none

    ! subroutine arguments
    type(receivers) :: recs
    type(receiver_mode):: recmode
    type(wavefront_container):: wafc
type(rat_conf):: conf

    ! local varaibles 
    integer :: i,j



    ! reset wavefronts storage
    do i=1,wafc%n
       if (wafc%waf(i)%n_nod>0) then

          deallocate(wafc%waf(i)%pos)
          deallocate(wafc%waf(i)%spf)
          deallocate(wafc%waf(i)%neig)
          deallocate(wafc%waf(i)%con)
          deallocate(wafc%waf(i)%stat)
          deallocate (wafc%waf(i)%index%ind)
          wafc%waf(i)%index%n=0
          wafc%waf(i)%n_nod=0

       end if
    end do

    ! reset receiver storage

if (conf%do_rays>0) then
    do i=1,recs%nn
       do j=1,recs%nar(i)

          recs%ariv(i,j)%path%n=0
          deallocate (recs%ariv(i,j)%path%pos)

       end do
    end do
end if
    recs%nar=0

    !recs%artime=0.0
    !recs%path=0.0
    !recs%paic=0
    !recs%arind=0
    !recs%arpid=0
    !    recs%rel=0.0





    ! setup the receiver
    if (conf%recmode==1) then
       recs%stat=0
       do i=1,recmode%act(conf%sid)%n
          recs%stat(recmode%act(conf%sid)%val(i))=1
       end do
    else 
       recs%stat=1
    end if




  end subroutine reset_wavefronts_and_arrivals

!-----------------------------------------------------------------------------------

  subroutine update_strip(head,vmod,conf)


    use my_types
    use my_functions
    use my_constants
    use m_solver

    implicit none

    ! subroutine arguments
    type(rat_conf)::conf
    type(node),pointer :: head
    type(scalar_field_2d) :: vmod

    !local variables
    type(node),pointer ::cur=>null()
    real(kind=dbl)::x(3),xout(3),xerr(3)
    real(kind=dbl):: eps,htry,hmin

    cur=>head
    do
       if (.not.associated(cur)) exit  ! if we have update the tail

       x(1)=cur%pos(1)
       x(2)=cur%pos(2)
       x(3)=cur%pos(3)

       select case(conf%wt%solver)
       case(1)
          call rk1_k(x,xout,conf%wt%dt,vmod)
       case(2)
          call rk4_k(x,xout,conf%wt%dt,vmod)
       case(3)
          call rkck5_k(x,xout,xerr,conf%wt%dt,vmod)
       case(4)
          htry=conf%wt%dt         ! we try the default time step
          eps=10**(-6)          ! desired accuracy
          hmin=conf%wt%dt/10.0    ! minimum tiem step size
          call odeint_k(x,xout,conf%wt%dt,eps,hmin,vmod)
       end select

       cur%pos(1)=xout(1)
       cur%pos(2)=xout(2)
       cur%pos(3)=xout(3)

       cur=>cur%next
    end do

  end subroutine update_strip

  !--------------------------------------------------------------------------------!

  subroutine update_spfa(head,tail)

    use my_types 
    use my_functions
    implicit none

    type(node),pointer :: head,tail
    type(node),pointer ::cur=>null()

    real(kind=dbl):: d1,d2,d



    if (head%con(2)==1) then
       d1=great_circle_distance(head%pos(1),head%next%pos(1),head%pos(2),head%next%pos(2))/2.0
    else
       d1=0
    end if

    cur=>head%next
    do
       if (.not.associated(cur%next)) exit  ! if we have update the tail
       if (cur%con(2)==1) then
          d2=great_circle_distance(cur%pos(1),cur%next%pos(1),cur%pos(2),cur%next%pos(2))/2.0
       else
          d2=0
       end if
       d=d1+d2
  
             cur%spf=cur%spf+(d-cur%dinr)
   
       d1=d2
       cur=>cur%next
    end do


! check head ray/node

  if (tail%con(2)==1) then
       d1=great_circle_distance(head%pos(1),tail%pos(1),head%pos(2),tail%pos(2))/2.0
    else
       d1=0
    end if

    if (head%con(2)==1) then
       d2=great_circle_distance(head%pos(1),head%next%pos(1),head%pos(2),head%next%pos(2))/2.0
    else
       d2=0
    end if

    d=d1+d2


       head%spf=head%spf+(d-head%dinr)
    


! chekc tail ray/node

    if (tail%con(2)==1) then
       d1=great_circle_distance(head%pos(1),tail%pos(1),head%pos(2),tail%pos(2))/2.0
    else
       d1=0
    end if

    if (tail%prev%con(2)==1) then
       d2=great_circle_distance(tail%pos(1),tail%prev%pos(1),tail%pos(2),tail%prev%pos(2))/2.0
    else
       d2=0
    end if
d=d1+d2

      tail%spf=tail%spf+(d-tail%dinr)
    




  end subroutine update_spfa


  !--------------------------------------------------------------------------------!


  subroutine update_dinr(head,tail)

    use my_types
    use my_functions
    implicit none

    type(node),pointer :: head,tail
    type(node),pointer ::cur=>null()

    real(kind=dbl):: d1,d2

    if (head%con(2)==1) then
       d1=great_circle_distance(head%pos(1),head%next%pos(1),head%pos(2),head%next%pos(2))/2.0
    else
       d1=0
    end if

    cur=>head%next
    do
       if (.not.associated(cur%next)) exit  ! if we have update the tail
       if (cur%con(2)==1) then
          d2=great_circle_distance(cur%pos(1),cur%next%pos(1),cur%pos(2),cur%next%pos(2))/2.0
       else 
          d2=0
       end if
       cur%dinr=d1+d2
       d1=d2
       cur=>cur%next
    end do

    if (tail%con(2)==1) then
       d1=great_circle_distance(head%pos(1),tail%pos(1),head%pos(2),tail%pos(2))/2.0
    else
       d1=0
    end if

    if (head%con(2)==1) then
       d2=great_circle_distance(head%pos(1),head%next%pos(1),head%pos(2),head%next%pos(2))/2.0
    else
       d2=0
    end if

    head%dinr=d1+d2



    if (tail%con(2)==1) then
       d1=great_circle_distance(head%pos(1),tail%pos(1),head%pos(2),tail%pos(2))/2.0
    else
       d1=0
    end if

    if (tail%prev%con(2)==1) then    
       d2=great_circle_distance(tail%pos(1),tail%prev%pos(1),tail%pos(2),tail%prev%pos(2))/2.0
    else
       d2=0
    end if
    tail%dinr=d1+d2


  end subroutine update_dinr

  !---------------------------------------------------------------------------!

  subroutine insert_nodes_linear(head,tail,params)

    ! If the pahse space distance between two points on the bicharacteristic strip becomes 
    ! bigger than a predefined threshold a point is interpolated using a linear itnerpolation

    use my_types
    use my_functions
    use my_constants
    use m_solver
    implicit none

    ! subroutine arguments
    type(rat_conf)           :: params
    type(node),pointer    :: head,tail,cur

    !local variables
    type(node),pointer::newnode
    real(kind=dbl)::x1,y1,z1,x2,y2,z2
    real(kind=dbl)::d

    !write(*,*) 'insert nodes'

    ! isnert point between points which are to far apart from each other
    cur=>head
    do
       if (.not.associated(cur%next)) exit

       x1=cur%pos(1)
       y1=cur%pos(2)
       z1=cur%pos(3)

       x2=cur%next%pos(1)
       y2=cur%next%pos(2)
       z2=cur%next%pos(3)

       d=sqrt(((x1-x2)*params%wt%scax)**2+((y1-y2)*params%wt%scay)**2+(z1-z2)**2)

       if ((cur%con(2)==1).and.(cur%next%con(1)==1)) then
          if (d>params%wt%pdma) then
             ! insert points
             allocate(newnode)
             newnode%pos(1)=(x1+x2)/2.0
             newnode%pos(2)=(y1+y2)/2.0
             newnode%pos(3)=(z1+z2)/2.0
             newnode%con(1)=1
             newnode%con(2)=1
             newnode%pid =params%wt%nfrid
             newnode%prev=>cur
             newnode%next=>cur%next
             cur%next%prev=>newnode
             cur%next=>newnode
             params%wt%n_nod=params%wt%n_nod+1
             params%wt%nfrid=params%wt%nfrid+1
             cur=>head

             newnode%spf=(cur%spf+cur%next%spf)/2.0

          end if
       end if
       cur=>cur%next
    end do


    ! insert point between head and tail if they are to far appart

    x1=head%pos(1)
    y1=head%pos(2)
    z1=head%pos(3)+2.0*pi
    x2=tail%pos(1)
    y2=tail%pos(2)
    z2=tail%pos(3)

    d=sqrt(((x1-x2)*params%wt%scax)**2+((y1-y2)*params%wt%scay)**2+(z1-z2)**2)

    if (head%con(1)==1.and.tail%con(2)==1) then
       if (d>params%wt%pdma) then
          allocate(tail%next)
          tail%next%prev=>tail
          tail => tail%next
          nullify(tail%next)          
          tail%pos(1)=(x1+x2)/2.0
          tail%pos(2)=(y1+y2)/2.0
          tail%pos(3)=(z1+z2)/2.0
          tail%con(1)=1
          tail%con(2)=1
          tail%pid=params%wt%nfrid 
          params%wt%n_nod=params%wt%n_nod+1
          params%wt%nfrid=params%wt%nfrid+1
          tail%spf=(head%spf+tail%spf)/2.0
       end if
    end if

  end subroutine insert_nodes_linear

  !-----------------------------------------------------------------------------!


  subroutine insert_nodes_splines(head,tail,params)

    ! If the pahse space distance between two points on the bicharacteristic strip becomes 
    ! bigger than a predefined threshold a point is interpolated using a linear itnerpolation

    use my_types
    use my_functions
    use my_constants
    use m_solver
    implicit none

    ! subroutine arguments
    type(rat_conf)           :: params
    type(node),pointer    :: head,tail,cur

    !local variables
    type(node),pointer::newnode
    real(kind=dbl)::x1,y1,z1,x2,y2,z2
    real(kind=dbl):: x0,y0,z0,x3,y3,z3
    real(kind=dbl)::d

    !write(*,*) 'insert nodes'

    ! isnert point between points which are to far apart from each other
    cur=>head
    do
       if (.not.associated(cur%next)) exit

       x1=cur%pos(1)
       y1=cur%pos(2)
       z1=cur%pos(3)

       x2=cur%next%pos(1)
       y2=cur%next%pos(2)
       z2=cur%next%pos(3)

       d=sqrt(((x1-x2)*params%wt%scax)**2+((y1-y2)*params%wt%scay)**2+(z1-z2)**2)

       if ((cur%con(2)==1).and.(cur%next%con(1)==1)) then
          if (d>params%wt%pdma) then
             ! insert points
             allocate(newnode)

             ! we do the lienar interpoaltion and check afterwards if a spline is possible.
             ! if we can use a splien we do it using a spline

             newnode%pos(1)=(x1+x2)/2.0
             newnode%pos(2)=(y1+y2)/2.0
             newnode%pos(3)=(z1+z2)/2.0

             if (associated(cur%prev).and.associated(cur%next%next)) then
                if (cur%con(1)==1.and.cur%next%con(2)==1) then

                   x0=cur%prev%pos(1)
                   y0=cur%prev%pos(2)
                   z0=cur%prev%pos(3)

                   x3=cur%next%next%pos(1)
                   y3=cur%next%next%pos(2)
                   z3=cur%next%next%pos(3)

                   newnode%pos(1)=1.0/16.0*(-1.0*x0+9.0*x1+9.0*x2-1.0*x3)
                   newnode%pos(2)=1.0/16.0*(-1.0*y0+9.0*y1+9.0*y2-1.0*y3)
                   newnode%pos(3)=1.0/16.0*(-1.0*z0+9.0*z1+9.0*z2-1.0*z3)

                end if
             end if


             newnode%con(1)=1
             newnode%con(2)=1
             newnode%pid =params%wt%nfrid
             newnode%prev=>cur
             newnode%next=>cur%next
             cur%next%prev=>newnode
             cur%next=>newnode
             params%wt%n_nod=params%wt%n_nod+1
             params%wt%nfrid=params%wt%nfrid+1
             cur=>head
          end if
       end if
       cur=>cur%next
    end do


    ! insert point between head and tail if they are to far appart

    x1=head%pos(1)
    y1=head%pos(2)
    z1=head%pos(3)+2.0*pi
    x2=tail%pos(1)
    y2=tail%pos(2)
    z2=tail%pos(3)

    d=sqrt(((x1-x2)*params%wt%scax)**2+((y1-y2)*params%wt%scay)**2+(z1-z2)**2)

    if (head%con(1)==1.and.tail%con(2)==1) then
       if (d>params%wt%pdma) then

          newnode%pos(1)=(x1+x2)/2.0
          newnode%pos(2)=(y1+y2)/2.0
          newnode%pos(3)=(z1+z2)/2.0

          if (associated(tail%prev).and.associated(head%next)) then
             if (tail%con(1)==1.and.head%con(2)==1) then

                x0=head%next%next%pos(1)
                y0=head%next%next%pos(2)
                z0=head%next%next%pos(3)+2.0*pi

                x3=tail%prev%pos(1)
                y3=tail%prev%pos(2)
                z3=tail%prev%pos(3)

                newnode%pos(1)=1.0/16.0*(-1.0*x0+9.0*x1+9.0*x2-1.0*x3)
                newnode%pos(2)=1.0/16.0*(-1.0*y0+9.0*y1+9.0*y2-1.0*y3)
                newnode%pos(3)=1.0/16.0*(-1.0*z0+9.0*z1+9.0*z2-1.0*z3)

             end if
          end if



          allocate(tail%next)
          tail%next%prev=>tail
          tail => tail%next
          nullify(tail%next)          
          tail%pos(1)=(x1+x2)/2.0
          tail%pos(2)=(y1+y2)/2.0
          tail%pos(3)=(z1+z2)/2.0
          tail%con(1)=1
          tail%con(2)=1
          params%wt%n_nod=params%wt%n_nod+1
          tail%pid=params%wt%nfrid
          params%wt%nfrid=params%wt%nfrid+1
       end if
    end if

  end subroutine insert_nodes_splines

  !-----------------------------------------------------------------------------!

  subroutine remove_nodes(head,tail,params)

    ! If two points on the bicharacteristic strip are closer together than a predefined
    !  threshold. The poitn between thes two point is removed. Poitns whcih ar eleaving 
    ! the spatial domain are also remcoed from the bicharacteristic strip

    use my_types
    use my_functions
    use my_constants

    use m_solver

    implicit none

    ! subroutine arguments
    type(rat_conf)                   :: params
    type(node),pointer                :: head,tail

    !local variables
    type(node),pointer :: oldnode=>null()
    type(node),pointer :: cur=>null()
    real(kind=dbl)     :: x1,y1,z1,x2,y2,z2,x,y,z,d
    integer :: switch


    !write(*,*) 'remove nodes'

    ! remove poitns outside of the computational domain
    cur=>head
    do
       if (.not.associated(cur)) exit

       switch=0

       x=cur%pos(1)
       y=cur%pos(2)
       z=cur%pos(3)
       if (((x<=params%xmin.or.&
            (x>=params%xmax)).or.&
            ((y<=params%ymin).or.&
            (y>=params%ymax))))  then
          switch=1
       end if

       if (associated(cur%prev)) then
          if (cur%con(1)==1) then
             x=cur%prev%pos(1)
             y=cur%prev%pos(2)
             z=cur%prev%pos(3)
             if (((x<=params%xmin.or.&
                  (x>=params%xmax)).or.&
                  ((y<=params%ymin).or.&
                  (y>=params%ymax))))  then
                switch=switch+1
             end if
          end if
       end if


       if (associated(cur%next)) then
          if (cur%con(2)==1) then
             x=cur%next%pos(1)
             y=cur%next%pos(2)
             z=cur%next%pos(3)

             if (((x<=params%xmin.or.&
                  (x>=params%xmax)).or.&
                  ((y<=params%ymin).or.&
                  (y>=params%ymax))))  then



                switch=switch+1
             end if
          end if
       end if



       if (switch>=2) then

          if (.not. associated(cur%next)) then
             ! remove tail
             oldnode => tail
             tail => tail%prev
             tail%con(1)=1
             tail%con(2)=0
             head%con(1)=0
             nullify(tail%next)
             deallocate(oldnode)
             cur => head
             params%wt%n_nod=params%wt%n_nod-1
          else if (.not. associated (cur%prev)) then
             ! remove head
             oldnode => head
             head=>head%next
             head%con(1)=0
             head%con(2)=1
             tail%con(2)=0
             nullify(head%prev)
             deallocate(oldnode)
             cur=>head
             params%wt%n_nod=params%wt%n_nod-1
          else
             ! remove in between
             oldnode => cur
             cur%prev%next => cur%next
             cur%next%prev => cur%prev
             cur%prev%con(2)=0
             cur%next%con(1)=0
             deallocate(oldnode)
             params%wt%n_nod=params%wt%n_nod-1
             cur=>head
          end if
          ! it makes no sense to evolve a bicharacteristic strip with three points on it
          if (params%wt%n_nod<=3) then
             return
          end if
       end if

       if (params%wt%n_nod<=3) then
          return
       end if

       cur=>cur%next
    end do

    ! remove point between two points which are to close together
    cur=>head%next
    do
       if (.not.associated(cur%next)) exit

       x1=cur%prev%pos(1)
       y1=cur%prev%pos(2)
       z1=cur%prev%pos(3)

       x2=cur%next%pos(1)
       y2=cur%next%pos(2)
       z2=cur%next%pos(3)

       d=sqrt(((x1-x2)*params%wt%scax)**2+((y1-y2)*params%wt%scay)**2+(z1-z2)**2)

       if (d<params%wt%pdmi) then
          ! remove point
          oldnode => cur
          cur%prev%next => cur%next
          cur%next%prev => cur%prev
          cur => head%next    
          deallocate(oldnode)
          params%wt%n_nod=params%wt%n_nod-1
       end if

       ! it makes no sense to evovle a bicharacteristic strip with three points on it
       if (params%wt%n_nod<=3) then
          return
       end if

       cur=>cur%next
    end do

  end subroutine remove_nodes

  !---------------------------------------------------------------------------------!

  subroutine store_strip(head,tail,wafc,conf)

    ! This subroutine stores tha bciharactersitic strip and the connectivity between
    ! the points on the strip in order to acces this information alter when doing the ray 
    ! tracing from the receiver back to the source for the different arrivals

    use my_types
    use my_functions

    use m_solver

    implicit none

    ! subroutine arguments
    type(rat_conf)            :: conf
    type(node),pointer    :: head,tail,cur
    !  type(wavefronts)           :: wafs
    type(wavefront_container):: wafc
    integer ::i,m


    wafc%waf(conf%wt%curit)%n_nod=conf%wt%n_nod
    allocate(wafc%waf(conf%wt%curit)%pos(wafc%waf(conf%wt%curit)%n_nod,3))
    allocate(wafc%waf(conf%wt%curit)%spf(wafc%waf(conf%wt%curit)%n_nod))
    allocate(wafc%waf(conf%wt%curit)%neig(wafc%waf(conf%wt%curit)%n_nod,2))
    allocate(wafc%waf(conf%wt%curit)%con(wafc%waf(conf%wt%curit)%n_nod,2))
    allocate(wafc%waf(conf%wt%curit)%stat(wafc%waf(conf%wt%curit)%n_nod))

    wafc%waf(conf%wt%curit)%index%n=conf%wt%n_nod
    allocate (wafc%waf(conf%wt%curit)%index%ind(wafc%waf(conf%wt%curit)%index%n,2))


    i=1
    cur=>head

    do

       wafc%waf(conf%wt%curit)%pos(i,1)=cur%pos(1)
       wafc%waf(conf%wt%curit)%pos(i,2)=cur%pos(2)
       wafc%waf(conf%wt%curit)%pos(i,3)=cur%pos(3)
       wafc%waf(conf%wt%curit)%spf(i)=cur%spf



       if ((associated(cur%prev)).and.(associated(cur%next))) then
          wafc%waf(conf%wt%curit)%neig(i,1)=cur%prev%pid
          wafc%waf(conf%wt%curit)%neig(i,2)=cur%next%pid
          wafc%waf(conf%wt%curit)%con(i,1)=cur%con(1)
          wafc%waf(conf%wt%curit)%con(i,2)= cur%con(2)
       else if ( .not. associated(cur%prev)) then
          wafc%waf(conf%wt%curit)%neig(i,1)=tail%pid
          wafc%waf(conf%wt%curit)%neig(i,2)=head%next%pid
          wafc%waf(conf%wt%curit)%con(i,1)=cur%con(1)
          wafc%waf(conf%wt%curit)%con(i,2)= cur%con(2)              
       else if (.not. associated(cur%next)) then       
          wafc%waf(conf%wt%curit)%neig(i,1)=tail%prev%pid
          wafc%waf(conf%wt%curit)%neig(i,2)=head%pid
          wafc%waf(conf%wt%curit)%con(i,1)=cur%con(1)
          wafc%waf(conf%wt%curit)%con(i,2)= cur%con(2)   
       end if

       wafc%waf(conf%wt%curit)%index%ind(i,:)=(/cur%pid,i/)
       !        if (conf%wt%curit==500) then
       !write(13,*) cur%pid,i
       !end if
       if (.not.associated(cur%next)) exit
       cur=>cur%next
       i=i+1


    end do




    call heap_sort_list_index(wafc%waf(conf%wt%curit)%index)
    wafc%waf(conf%wt%curit)%headpid=head%pid







  end subroutine store_strip





end module m_wavefronts



