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
! * m_arrivals *
!
! This module contains subroutines used to extract and compute arrival information
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
!   You may contact the author by:0
!       e-mail:  juerg@rses.anu.edu.au
!

module m_arrivals

  implicit none

contains

  !-------------------------------------------------------------------------------------!

  subroutine receiver_processing(recs,wafc,conf)
    ! computing the exact arrival time. The recorded squard distanc eto the closest
    ! point on thw wavefrotn is used to trigge the more advanced implementation which 
    ! should recognize the passige of the wavefornt segmetns in the same time step.
    ! the code

    !    - -------------
    !     |          |
    !   k2 h2  r   k1 h1
    !     |          |
    !    ---------------
    !

    use my_types
    use my_constants
    use my_functions
    implicit none

    !subroutine arguments
    type(receivers)            :: recs
    type(wavefront_container)  :: wafc
    type(rat_conf)             :: conf

    ! local variables
    integer :: i,j,m,pol,n,mm
    integer :: k1,k2,k1_
    integer  :: h1,h2,h1_
    logical :: bol1,bol2
  !  real(kind=dbl):: a,b,c
    real(kind=dbl):: x(4),y(4),z(4),xr,yr,xc1,yc1,xc2,yc2,zc1,zc2,d0,d1,d2,d3,d4,d5
    real(kind=dbl):: s(4),sc1,sc2
    real(kind=dbl) :: artime(recs%mar+3),arazi(recs%mar+3),arspf(recs%mar+3)
    integer :: arpid(recs%mar+3),arind(recs%mar+3),ind



    ! a potential passage of the wavefront at the receiver testing the segments

    do i=1,recs%nn                        ! loop over all receivers
       recs%nar(i)=0
   if (recs%stat(i)==1) then
      !write(*,*) recs%pos(i,:)

       artime=conf%wt%maxit*4*conf%wt%dt  ! initialise arrival time vector
       k1_=-1  

       j=2   ! We do not use the wavefront at the source because nobody healthy 
       ! in his mind puts a receiver on top of a source - except es


       ! there is not much point in using wavefronts which consists of less 
       ! than 10 nodes  for searching an arrival along them...

       do while(j<conf%wt%maxit.and.wafc%waf(j)%n_nod>10)

          ! we are searching for the recevier only in front of the wavefront and 
          ! not behind the wavefront as this search is done in the previous step.

          k1=wafc%waf(j)%headpid

          m=1
          ind=0 !counter for the arrivals of this wavefront

          ! initailising the arrival time vector for the current time step and receiver

          artime=conf%wt%maxit*2*conf%wt%dt
          ! set the index of the previously used right ray to a non exsistent ray
          k1_=-1   

          do while (m<=wafc%waf(j)%n_nod)

             call locate_in_list_index(wafc%waf(j)%index,k1,h1)
             k2=wafc%waf(j)%neig(wafc%waf(j)%index%ind(h1,2),2)
             call locate_in_list_index(wafc%waf(j)%index,k2,h2)


             mm=0
             pol=1  ! we assume that we can construct a polygon

             ! We will now perform diferent checks if we really can
             ! construct a polygon

             ! polygon construction not possible, there is no connection 
             ! between k1 and k2
             if (wafc%waf(j)%con(wafc%waf(j)%index%ind(h1,2),2)==0) then
                pol=0
             end if

             ! making sure that a polygon can be constructed using the current
             ! ray segments. Whcih means k1 and ks have to presented also at time i+1
             ! otherwise we take the left respecteviely right neighbour of k1 
             ! respectively k2 if we can take the neighbours due to a lack of 
             ! connectivity the plygon construction is not possible and hence 
             ! pol=0 and the process is stopped

             bol1=is_in_list_index(wafc%waf(j)%index,k1)
             bol2=is_in_list_index(wafc%waf(j+1)%index,k1)

             do while (.not.(bol1).or..not.(bol2))
                if (wafc%waf(j)%con(wafc%waf(j)%index%ind(h1,2),1)==1) then
                   k1=wafc%waf(j)%neig(wafc%waf(j)%index%ind(h1,2),1)
                   call locate_in_list_index(wafc%waf(j)%index,k1,h1)
                else 
                   pol=0
                   exit
                end if
                bol1=is_in_list_index(wafc%waf(j)%index,k1)
                bol2=is_in_list_index(wafc%waf(j+1)%index,k1)
             end do

             bol1=is_in_list_index(wafc%waf(j)%index,k2)
             bol2=is_in_list_index(wafc%waf(j+1)%index,k2)

             do while (.not.(bol1).or..not.(bol2))
                if (wafc%waf(j)%con(wafc%waf(j)%index%ind(h2,2),2)==1) then
                   k2=wafc%waf(j)%neig(wafc%waf(j)%index%ind(h2,2),2)
                   call locate_in_list_index(wafc%waf(j)%index,k2,h2)
                   mm=mm+1
                else 
                   pol=0
                   exit
                end if
                bol1=is_in_list_index(wafc%waf(j)%index,k2)
                bol2=is_in_list_index(wafc%waf(j+1)%index,k2)
             end do

             ! mm is used for measuring the motion of the ray2 side of the polygon
             ! when constructing one.
             if (mm>0) then
                m=m+mm
             end if

             ! k1 is the same as the previous k1 hence the construction of a 
             ! polygon was impossible
             if (k1==k1_) then
                pol=0
             end if

             ! check if we can build the polygon based on the previous tests
             if (pol==1) then

                ! extracting the corners of the current polygon
                call locate_in_list_index(wafc%waf(j)%index,k1,h1)
                call locate_in_list_index(wafc%waf(j)%index,k2,h2)

                x(1)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h1,2),1)
                x(2)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h2,2),1)
                y(1)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h1,2),2)
                y(2)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h2,2),2)


                call locate_in_list_index(wafc%waf(j+1)%index,k1,h1)
                call locate_in_list_index(wafc%waf(j+1)%index,k2,h2)

                x(3)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h2,2),1)
                x(4)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h1,2),1)
                y(3)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h2,2),2)
                y(4)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h1,2),2)


!!$                call locate_in_list_index(wafc%waf(j)%index,k1,h1)
!!$                call locate_in_list_index(wafc%waf(j)%index,k2,h2)
!!$
!!$                y(1)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h1,2),2)
!!$                y(2)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h2,2),2)
!!$
!!$                call locate_in_list_index(wafc%waf(j+1)%index,k1,h1)
!!$                call locate_in_list_index(wafc%waf(j+1)%index,k2,h2)
!!$
!!$                y(3)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h2,2),2)
!!$                y(4)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h1,2),2)


                ! extracting the receiver position
                xr=recs%pos(i,1)
                yr=recs%pos(i,2)

                ! check if the receiver lies in the current polygon
                if (in_polygon(x,y,xr,yr)>=0) then


                   ! the receiver is in the current polygon computing the exact arrival time
                   ! if the maximum number of arrivals to be detected has been reached we
                   ! abort
                   if (recs%nar(i)>=recs%mar) then
                      exit
                   end if





                   call closest_point_on_line(xr,yr,x(1),y(1),x(2),y(2),xc1,yc1)
                   call closest_point_on_line(xr,yr,x(3),y(3),x(4),y(4),xc2,yc2)


                   d0=great_circle_distance(xc2,xc1,yc2,yc1)
                   d1=great_circle_distance(xr,xc1,yr,yc1)

                   ind=ind+1

                   artime(ind)=(j-1)*conf%wt%dt+(conf%wt%dt/d0*d1)

                   arind(ind)=j

                   d1=sqrt((x(1)-xr)**2+(y(1)-yr)**2)
                   d2=sqrt((x(2)-xr)**2+(y(2)-yr)**2)

                   if (minloc((/d1,d2/),1)==1) then
                      arpid(ind)=k1
                   else
                      arpid(ind)=k2
                   end if

                   call locate_in_list_index(wafc%waf(j)%index,k1,h1)
                   call locate_in_list_index(wafc%waf(j)%index,k2,h2)

                   z(1)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h1,2),3)
                   z(2)=wafc%waf(j)%pos(wafc%waf(j)%index%ind(h2,2),3)

                   call locate_in_list_index(wafc%waf(j+1)%index,k1,h1)
                   call locate_in_list_index(wafc%waf(j+1)%index,k2,h2)

                   z(3)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h2,2),3)
                   z(4)=wafc%waf(j+1)%pos(wafc%waf(j+1)%index%ind(h1,2),3)
                   ! azimuth determination

                   d2=great_circle_distance(x(1),x(2),y(1),y(2))
                   d3=great_circle_distance(x(1),xc1,y(1),yc1)

                   d4=great_circle_distance(x(3),x(4),y(3),y(4))
                   d5=great_circle_distance(x(3),xc2,y(3),yc2)

                   zc1=z(1)+(z(2)-z(1))*d3/d2     
                   zc2=z(3)+(z(4)-z(3))*d5/d4


                   arazi(ind)=(zc1+(zc2-zc1)*d1/d0)*radian
                   if (arazi(ind)>0.0) then
arazi(ind)=abs(180-arazi(ind)) 
else
arazi(ind)=abs(-180-arazi(ind))
end if
                   if (conf%wt%mode==2) then


                      call locate_in_list_index(wafc%waf(j)%index,k1,h1)
                      call locate_in_list_index(wafc%waf(j)%index,k2,h2)

                      s(1)=wafc%waf(j)%spf(wafc%waf(j)%index%ind(h1,2))
                      s(2)=wafc%waf(j)%spf(wafc%waf(j)%index%ind(h2,2))


                      call locate_in_list_index(wafc%waf(j+1)%index,k1,h1)
                      call locate_in_list_index(wafc%waf(j+1)%index,k2,h2)

                      s(3)=wafc%waf(j+1)%spf(wafc%waf(j+1)%index%ind(h2,2))
                      s(4)=wafc%waf(j+1)%spf(wafc%waf(j+1)%index%ind(h1,2))

                      sc1=s(1)+(s(2)-s(1))*d3/d2     
                      sc2=s(3)+(s(4)-s(3))*d5/d4



                      arspf(ind)=sc1+(sc2-sc1)*d1/d0
                     
!write(*,*) (arspf(ind))

                   end if


                end if

             end if ! in_polygon

             ! store the old k1 value
             k1_=k1
             !update the k value
             k1=k2
             m=m+1
          end do

          ! all arrivals of the current wavefront have been detected
          ! we need to sort them and store them with increasing arrival time



          do m=1,ind
             if (recs%nar(i)>=recs%mar) then
                exit
             end if
             recs%nar(i)=recs%nar(i)+1                      
             n=minloc(artime(1:ind),1)
             recs%ariv(i,recs%nar(i))%time=artime(n)
             recs%ariv(i,recs%nar(i))%iter=arind(n)
             recs%ariv(i,recs%nar(i))%raid=arpid(n)
             recs%ariv(i,recs%nar(i))%azi=arazi(n)
             recs%ariv(i,recs%nar(i))%spf=arspf(n)
             artime(n)=artime(n)*2.0


          end do
          j=j+1
       end do  ! end loop iterations
end if ! recs%stat(i)==1
    end do  ! end loop receivers



  end subroutine receiver_processing

  !-----------------------------------------------------------------------------------!

  subroutine extract_ray_paths_for_point_source(recs,wafc,conf)

    ! this subroutine is a wrap arround for the subroutine cosntruct_ray_path which extracts
    ! the ray path for the given arrival j at a given receiver i.
    ! The ray tracing is done for all arrivals ans all receivers the results are afterwards
    ! stored in recs

    use my_types
    implicit none

    ! subroutine arguments
    type(receivers)  :: recs
    type(wavefront_container):: wafc
    type(rat_conf)      :: conf

    ! local varaibles
    integer :: i,j,k
    real(kind=dbl):: toa,ain
    real(kind=dbl)::x(conf%wt%maxit),y(conf%wt%maxit)

    do i=1,recs%nn
       do j=1,recs%nar(i)
          !construct the ray path for the given receiver and arrival
          call construct_ray_path_for_point_source(i,j,wafc,recs,conf,x,y,toa,ain)

          recs%ariv(i,j)%path%n=recs%ariv(i,j)%iter+1
          allocate (recs%ariv(i,j)%path%pos(recs%ariv(i,j)%path%n,2))

          recs%ariv(i,j)%path%pos(:,1)=x(1:recs%ariv(i,j)%path%n)
          recs%ariv(i,j)%path%pos(:,2)=y(1:recs%ariv(i,j)%path%n)


       end do
    end do
  end subroutine extract_ray_paths_for_point_source

  !--------------------------------------------------------------------------------------!

  subroutine construct_ray_path_for_point_source (ii,jj,wafc,recs,conf,xx,yy,toa,ain)

    ! this subroutine  extracts  the ray path for the given arrival jj at a given receiver i.

    use my_types
    use my_functions
    implicit none

    ! subroutine arguments
    type(receivers),intent(in)     :: recs
    type(wavefront_container)      :: wafc
    type(rat_conf),intent(in)      :: conf
    integer,intent(in)             :: ii,jj
    real(kind=dbl),intent(out):: xx(conf%wt%maxit),yy(conf%wt%maxit),toa,ain

    ! local variables
    real(kind=dbl) :: x,y,xr,yr,x1,y1,z1,x2,y2,z2,x2a,y2a,z2a,x2b,y2b,z2b,x1a,y1a,x3,y3
    real(kind=dbl):: d1,d2,d2a,d2b,dr,d3
    integer ::ray1,ray2a,ray2b,ray2,ray1a,ray3
    real(kind=dbl):: artime
    integer :: arind,i,iii
    integer :: h1,h2,h1a,h2a,h3
    logical:: bol1,bol2,bol3
    ! extract the receiver position
    xr=recs%pos(ii,1)
    yr=recs%pos(ii,2)

    ! time for the j-th arrival and pid
    ! iteration arrival time and interpolated arrival time
    arind=recs%ariv(ii,jj)%iter
    artime=recs%ariv(ii,jj)%time

    ! determine the ray responsibel for the arrival
    ray1=recs%ariv(ii,jj)%raid

    ! determine the neighbouring rays of ray1
    call locate_in_list_index(wafc%waf(arind)%index,ray1,h1)
    ray2a=wafc%waf(arind)%neig(wafc%waf(arind)%index%ind(h1,2),1)
    ray2b=wafc%waf(arind)%neig(wafc%waf(arind)%index%ind(h1,2),2)

    ! check if back tracking is possible

    if ((ray2a==0).or.(ray2b)==0) then
       write(*,*) 'warning neighbouring ray is missing'
       write(*,*) 'ray path extraction not possible'
       write(*,'(a12,i3)') ' source id', conf%sid
       write(*,'(a12,i3)') ' receiver id',ii
       write(*,'(a12,i3)') ' arrival id',jj
       return
    end if

    ! determine between ray1 and which of its two neighbours the receiver is positioned
    ! compute the position of the rays for the time artime

    call ray_interp(ray1,artime,arind,wafc,conf,x1,y1,z1)
    call ray_interp(ray2a,artime,arind,wafc,conf,x2a,y2a,z2a)
    call ray_interp(ray2b,artime,arind,wafc,conf,x2b,y2b,z2b)


    ! the back tracking is done in normal space
    ! compute the distance from the receiver to the rays 
    ! the two closest rays are the neighbours

    d1=sqrt((xr-x1)**2+(yr-y1)**2)
    d2a=sqrt((xr-x2a)**2+(yr-y2a)**2)
    d2b=sqrt((xr-x2b)**2+(yr-y2b)**2)

    ! check which ray is closer ray 2a or ray2b
    ! we want to have the ray between ray1 and ray2 and not between ray2 and ray1
    ! as this makes things simpler to implement afterwards. however require shere the
    ! copy of ray2 to ray1 and vice versa if the ray is between ray2 and ray1
    if (d2a<d2b) then
       x2=x1; y2=y1; z2=z1
       x1=x2a;y1=y2a;z1=z2a
       d2=d1
       d1=d2a
       ray2=ray1
       ray1=ray2a
    else ! d2a>=d2b 
       ray2=ray2b
       x2=x2b;y2=y2b;z2=z2b
       d2=d2b  
    end if

    call locate_in_list_index(wafc%waf(arind)%index,ray2,h2)


    ! chose closest ray semgent sif te inital ray segments can not be used for back tracking

    bol1=is_in_list_index(wafc%waf(arind-1)%index,ray1)
    bol2=is_in_list_index(wafc%waf(arind)%index,ray1)
    bol3=is_in_list_index(wafc%waf(arind+1)%index,ray1)

    do while (.not.(bol1).or..not.(bol2).or..not.(bol3))
       call locate_in_list_index(wafc%waf(arind)%index,ray1,h1)
       ray1=wafc%waf(arind)%neig(wafc%waf(arind)%index%ind(h1,2),1)

       x1=wafc%waf(arind)%pos(h1,1)
       y1=wafc%waf(arind)%pos(h1,2)


       d1=sqrt((xr-x1)**2+(yr-y1)**2)
       d2=sqrt((xr-x2)**2+(yr-y2)**2)
       dr=d1/(d1+d2)

       bol1=is_in_list_index(wafc%waf(arind-1)%index,ray1)
       bol2=is_in_list_index(wafc%waf(arind)%index,ray1)
       bol3=is_in_list_index(wafc%waf(arind+1)%index,ray1)

    end do

    bol1=is_in_list_index(wafc%waf(arind-1)%index,ray2)
    bol2=is_in_list_index(wafc%waf(arind)%index,ray2)
    bol3=is_in_list_index(wafc%waf(arind+1)%index,ray2)

    do while (.not.(bol1).or..not.(bol2).or..not.(bol3))
       call locate_in_list_index(wafc%waf(arind)%index,ray2,h2)
       ray2=wafc%waf(arind)%neig(wafc%waf(arind)%index%ind(h2,2),2)

       x2=wafc%waf(arind)%pos(h2,1)
       y2=wafc%waf(arind)%pos(h2,2)  

       d1=sqrt((xr-x1)**2+(yr-y1)**2)
       d2=sqrt((xr-x2)**2+(yr-y2)**2)
       dr=d1/(d1+d2)

       bol1=is_in_list_index(wafc%waf(arind-1)%index,ray2)
       bol2=is_in_list_index(wafc%waf(arind)%index,ray2)
       bol3=is_in_list_index(wafc%waf(arind+1)%index,ray2)

    end do

    ! store the start point which is the receiver
    ! set the index for the arrival times to the first iteration in the time 
    ! direction to the source.

    if (artime<=(arind-1)*conf%wt%dt) then
       arind=arind-1
    else
       arind=arind
    end if


    ! start point for the ray is the receiver
    xx(arind+1)=xr;yy(arind+1)=yr

    !distance ratio between d1 and the distance betwee ray1 and ray2
    dr=d1/(d1+d2)

    call locate_in_list_index(wafc%waf(arind)%index,ray1,h1)
    call locate_in_list_index(wafc%waf(arind)%index,ray2,h2)

    ain= wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h1,2),3)+&
         dr*(wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h2,2),3)-&
         wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h1,2),3))

    do i=arind,1,-1

       call locate_in_list_index(wafc%waf(i)%index,ray1,h1)
       call locate_in_list_index(wafc%waf(i)%index,ray2,h2)


       ! verify if the rays are still there or i rays have been interpolated
       ! When walking abck from the reciever the interpoaltion of rays is unlikely however
       ! certain rays being removed is not unusual. Hence we first check for the removal of ray
       ! check if one of the neighbouring rays has been removed

       bol1=is_in_list_index(wafc%waf(i)%index,ray1)
       bol2=is_in_list_index(wafc%waf(i)%index,ray2)

       do while (.not.(bol1).or..not.(bol2))
          !      pause
          if (bol1.and.bol2) then
             ! do nothing both rays are still present
             ! write(*,*) 'a'
          else if (.not.(bol1).and.bol2) then

             ray1=wafc%waf(i)%neig(wafc%waf(i)%index%ind(h2,2),1)

             call locate_in_list_index(wafc%waf(i+1)%index,ray1,h1)
             x1=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1,2),1)
             y1=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1,2),2)

             call locate_in_list_index(wafc%waf(i+1)%index,ray2,h2)
             x2=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2,2),1)
             y2=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2,2),2)


             x=xx(i+1)
             y=yy(i+1)

             d1=sqrt((x-x1)**2+(y-y1)**2)
             d2=sqrt((x-x2)**2+(y-y2)**2)
             dr=d1/(d1+d2)

          else if (bol1.and..not.(bol2)) then

             ray2=wafc%waf(i)%neig(wafc%waf(i)%index%ind(h1,2),2)

             call locate_in_list_index(wafc%waf(i+1)%index,ray1,h1)
             x1=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1,2),1)
             y1=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1,2),2)

             call locate_in_list_index(wafc%waf(i+1)%index,ray2,h2)
             x2=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2,2),1)
             y2=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2,2),2)

             x=xx(i+1)
             y=yy(i+1)

             d1=sqrt((x-x1)**2+(y-y1)**2)
             d2=sqrt((x-x2)**2+(y-y2)**2)
             dr=d1/(d1+d2)


          else if (.not.(bol1).and..not.(bol2)) then 

             call locate_in_list_index(wafc%waf(i+1)%index,ray1,h1)
             call locate_in_list_index(wafc%waf(i+1)%index,ray2,h2)

             ray1=wafc%waf(i+1)%neig(wafc%waf(i+1)%index%ind(h1,2),2)
             ray2=wafc%waf(i+1)%neig(wafc%waf(i+1)%index%ind(h2,2),2)

             ray1a=ray1
             call locate_in_list_index(wafc%waf(i+1)%index,ray1a,h1a)

             ray1=wafc%waf(i)%neig(wafc%waf(i)%index%ind(h2,2),1)
             call locate_in_list_index(wafc%waf(i)%index,ray1,h1)

             x1a=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1a,2),1)
             y1a=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1a,2),2)

             x1=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1,2),1)
             y1=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h1,2),2)


             ray2a=ray2
             call locate_in_list_index(wafc%waf(i+1)%index,ray2a,h2a)

             ray2=wafc%waf(i)%neig(wafc%waf(i)%index%ind(h1,2),2)
             call locate_in_list_index(wafc%waf(i)%index,ray2,h2)


             x2a=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2a,2),1)
             y2a=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2a,2),2)

             x2=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2,2),1)
             y2=wafc%waf(i+1)%pos(wafc%waf(i+1)%index%ind(h2,2),2)


             x=xx(i+1)
             y=yy(i+1)

             d1=sqrt((x-x1)**2+(y-y1)**2)
             d2=sqrt((x-x2)**2+(y-y2)**2)
             dr=d1/(d1+d2)
             write(*,*) 'd'
          end if

          call locate_in_list_index(wafc%waf(i)%index,ray1,h1)
          call locate_in_list_index(wafc%waf(i)%index,ray2,h2)

          bol1=is_in_list_index(wafc%waf(i)%index,ray1)
          bol2=is_in_list_index(wafc%waf(i)%index,ray2)

       end do


       ! we need to test if thier is a new ray appearing it is iunlikely as we are 
       ! walking form the receiver to the source however you never know


       do while (wafc%waf(i)%neig(wafc%waf(i)%index%ind(h1,2),2)/= wafc%waf(i)%index%ind(h2,1))


          ! There is a ray appearing between ray1 and ray2
          ! we need to determine whether or not it is left or right of the current 
          ! ray point based on  ray1 and ray2. It is assumed that ray1,  ray2 and ray3
          ! lie on a line. This is nor entirely correct, however a good assumption 


          x1=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h1,2),1)
          y1=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h1,2),2)

          x2=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h2,2),1)
          y2=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h2,2),2)


          x=x1+dr*(x2-x1)
          y=y1+dr*(y2-y1)

          ray3=wafc%waf(i)%neig(wafc%waf(i)%index%ind(h1,2),2)
          call locate_in_list_index(wafc%waf(i)%index,ray3,h3)

          x3=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h3,2),1)
          y3=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h3,2),2)   


          d1=sqrt((x-x1)**2+(y-y1)**2)
          d3=sqrt((x3-x1)**2+(y3-y1)**2)


          if (d3>d1) then
             ray2=ray3
          else
             ray1=ray3
          end if

          call locate_in_list_index(wafc%waf(i)%index,ray1,h1)
          call locate_in_list_index(wafc%waf(i)%index,ray2,h2)


          x1=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h1,2),1)
          y1=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h1,2),2)


          x2=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h2,2),1)
          y2=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h2,2),2)

          d1=sqrt((x-x1)**2+(y-y1)**2)
          d2=sqrt((x-x2)**2+(y-y2)**2)
          dr=d1/(d1+d2)



       end do

       x1=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h1,2),1)
       y1=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h1,2),2)

       x2=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h2,2),1)
       y2=wafc%waf(i)%pos(wafc%waf(i)%index%ind(h2,2),2)

       xx(i)=x1+dr*(x2-x1)
       yy(i)=y1+dr*(y2-y1)

    end do

    toa= wafc%waf(1)%pos(wafc%waf(1)%index%ind(h1,2),3)+&
         dr*(wafc%waf(1)%pos(wafc%waf(1)%index%ind(h2,2),3)-&
         wafc%waf(1)%pos(wafc%waf(1)%index%ind(h1,2),3))

  end subroutine construct_ray_path_for_point_source

  !-----------------------------------------------------------------------------------!

  subroutine ray_interp(pid,artime,arind,wafc,conf,xx,yy,zz)

    ! This subroutine linear interpolates the position of point of the bicharacteristic 
    ! strip for at itme betwwen tow discret tiem steps.

    use my_types
    use my_functions
    implicit none

    ! subroutine arguments
    integer    :: pid,arind
    real(kind=dbl) :: artime
    type(wavefront_container):: wafc
    real(kind=dbl):: xx,yy,zz
    type(rat_conf):: conf
    ! local variables

    real(kind=dbl) :: x1,y1,z1,x2,y2,z2
    real(kind=dbl)::t1,t2
    integer :: h1,h2

    ! determine if the interpolated arrival time is before or after discrete 
    ! time step arrival time

    if (artime<=(arind-1)*conf%wt%dt) then


       call locate_in_list_index(wafc%waf(arind-1)%index,pid,h1)
       call locate_in_list_index(wafc%waf(arind)%index,pid,h2)


       x1=wafc%waf(arind-1)%pos(wafc%waf(arind-1)%index%ind(h1,2),1)
       y1=wafc%waf(arind-1)%pos(wafc%waf(arind-1)%index%ind(h1,2),2)
       z1=wafc%waf(arind-1)%pos(wafc%waf(arind-1)%index%ind(h1,2),3)
       t1=(arind-1-1)*conf%wt%dt

       x2=wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h2,2),1)
       y2=wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h2,2),2)
       z2=wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h2,2),3)
       t2=(arind-1)*conf%wt%dt

    else 

       call locate_in_list_index(wafc%waf(arind)%index,pid,h1)
       call locate_in_list_index(wafc%waf(arind+1)%index,pid,h2)


       x1=wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h1,2),1)
       y1=wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h1,2),2)
       z1=wafc%waf(arind)%pos(wafc%waf(arind)%index%ind(h1,2),3)
       t1=(arind-1)*conf%wt%dt


       x2=wafc%waf(arind+1)%pos(wafc%waf(arind+1)%index%ind(h2,2),1)
       y2=wafc%waf(arind+1)%pos(wafc%waf(arind+1)%index%ind(h2,2),2)
       z2=wafc%waf(arind+1)%pos(wafc%waf(arind+1)%index%ind(h2,2),3)
       t2=(arind+1-1)*conf%wt%dt

    end if

    call raypath_linear_interp_3d(t1,x1,y1,z1,t2,x2,y2,z2,artime,xx,yy,zz)


  end subroutine ray_interp

  !--------------------------------------------------------------------------------!

  subroutine analyze_ray_paths_for_point_source(recs,conf)

    ! this subroutine is a wrap arround for the subroutine cosntruct_ray_path which extracts
    ! the ray path for the given arrival j at a given receiver i.
    ! The ray tracing is done for all arrivals ans all receivers the results are afterwards
    ! stored in recs

    use my_types
    implicit none

    ! subroutine arguments
    type(receivers) :: recs
    type(rat_conf) :: conf

    ! local varaibles
    integer :: i,j

    do i=1,recs%nn
       do j=1,recs%nar(i)
          if ((minval(recs%ariv(i,j)%path%pos(1:recs%ariv(i,j)%path%n,1))<= conf%xmin .or.&
               maxval(recs%ariv(i,j)%path%pos(1:recs%ariv(i,j)%path%n,1))>=conf%xmax).or.&
               (minval(recs%ariv(i,j)%path%pos(1:recs%ariv(i,j)%path%n,2))<= conf%ymin .or.&
               maxval(recs%ariv(i,j)%path%pos(1:recs%ariv(i,j)%path%n,2))>=conf%ymax)) then
             recs%ariv(i,j)%paic=-1 
          else 
             recs%ariv(i,j)%paic=1
             conf%tonar(j)=conf%tonar(j)+1
          end if
       end do
    end do

  end subroutine analyze_ray_paths_for_point_source


  !-------------------------------------------------

  subroutine setup_frechet_header(freha,vmod)
    ! Setups the header for a frechet matrix
    use my_types
    implicit none

    ! subroutine arguments
    type(frechet_matrix_header):: freha
    type(scalar_field_2d):: vmod

    ! local variables
    integer :: i,j, nn
    real(kind=dbl),pointer:: frerow(:)

    ! number of nodes
    freha%n2=(vmod%nx+vmod%cn*2)*(vmod%ny+vmod%cn*2)

    ! setup the number of rays
    freha%n1=0

    ! set the number of records
    freha%nre=0

    inquire(iolength=freha%recole) idum,idum,rdum

  end subroutine setup_frechet_header
  !--------------------------------------------------------------------------------!
  subroutine setup_frechet_matrix_files(vmod,conf)

    use my_types
    use my_functions
    use m_inout
    implicit none

    ! subroutine arguments
    type(Frechet_matrix_header):: freha
    type(rat_conf)             :: conf
    type(scalar_field_2d):: vmod
    ! local variables

    call  setup_frechet_header(freha,vmod)

    call write_frechet_header(freha,conf%ofn_frechet_hdr)

    call create_frechet_files(freha,conf%ofn_frechet_mat,conf%ofn_frechet_rai)

  end subroutine setup_frechet_matrix_files

  !----------------------------------------------------------------------------!

  subroutine save_frechet_matrix(recs,vmod,conf)

    use my_types
    use m_inout
    implicit none

    ! subroutine arguments
    type(rat_conf)             :: conf
    type(scalar_field_2d)        :: vmod
    type(receivers)            :: recs

    ! local variables
    type(frechet_matrix_header):: freha
    integer :: i,j,h,k
    integer :: iray
    type(frechet_matrix_row) :: frero
    integer :: n_frerow



    ! load the header
    call read_frechet_header(freha,conf%ofn_frechet_hdr)

    ! setup the current row for the frechet matrix...
    allocate(frero%val(freha%n2))
    allocate(frero%nze(freha%n2))

    ! we build the frechet matrix row by row that means raypath by rapyath...

    !loop ove rall the raypaths...

    do i=1,recs%nn
       do j=1,recs%nar(i)
          freha%n1=freha%n1+1
          frero%val=0.0_dbl;frero%nze=0
          call calc_frechet_row(vmod,recs%ariv(i,j)%path,freha,frero)
          call append_frechet_row(frero,freha,conf%ofn_frechet_mat,conf%ofn_frechet_rai)
       end do
    end do


    deallocate(frero%val,frero%nze)

    call write_frechet_header(freha,conf%ofn_frechet_hdr)

  end subroutine save_frechet_matrix



  !----------------------------------------------------------------------------!

  subroutine calc_frechet_row(vmod,path,freha,frero)

    use my_types
    implicit none

    ! subroutine arguments
    type(scalar_field_2d)     :: vmod
    type(raypath) :: path
    type(frechet_matrix_row):: frero
    type(frechet_matrix_header):: freha

    ! local variables
    integer :: k


    ! looping allong the raypath

!write(*,*) frero%val

!pause
    do k=2,path%n
       call calc_frechet_pathseg(path%pos(k-1,1),&
            path%pos(k-1,2),path%pos(k,1),path%pos(k,2),vmod,frero)    
    end do

!write(*,*) frero%nze
!write(*,*) frero%val
!pause
    !   close(13)

  end subroutine calc_frechet_row
  !--------------------------------------------------------------------------------!

  subroutine calc_frechet_pathseg(x1,y1,x2,y2,vmod,frero)

    use my_types
    use my_functions
    implicit none

    !subroutine arguments
    real(kind=dbl),intent(in):: x1,y1,x2,y2
    type(scalar_field_2d),intent(in):: vmod
    type(frechet_matrix_row)::frero
    ! local variables
    integer :: p1,q1,p2,q2
    real(kind=dbl):: ras(3,4)
    integer :: pq(3,2)
    integer :: nras  ! number of ray segments
    integer :: nodid
    real(kind=dbl) :: xa,ya,xb,yb,as(4),at(4),bs(4),bt(4),va,vb

    integer k,i,j,ii,jj
    real(kind=dbl):: ra,rb,dl
    ! determine for the point pair through which cell the connection is going

    ! localize the cell for x1 y1 
    p1=int((x1-vmod%x0)/vmod%dx)+1+vmod%cn
    q1=int((y1-vmod%y0)/vmod%dy)+1+vmod%cn

    ! localize the cell for x2 y2 
    p2=int((x2-vmod%x0)/vmod%dx)+1+vmod%cn
    q2=int((y2-vmod%y0)/vmod%dy)+1+vmod%cn

    ! check if they are in the same cell
    ! a ray segmetn can be extended over three cells
    ! check the segmentation of the ray
    if (p1==p2.and.q1==q2) then
       ! we have one ray segment
       ras(1,:)=(/x1,y1,x2,y2/)
       pq(1,:)=(/p1,q1/)
       nras=1
    end if

    if(p1/=p2.and.q1==q2) then
       ! we have two ray segments
       ! we go in the y direction
       xa=(max(p1,p2)-1-vmod%cn)*vmod%dx+vmod%x0
       ya=linear_interp_1d(x1,y1,x2,y2,xa)
       ras(1,:)=(/x1,y1,xa,ya/)
       ras(2,:)=(/xa,ya,x2,y2/)
       pq(1,:)=(/p1,q1/)
       pq(2,:)=(/p2,q2/)
       nras=2
    end if

    if(p1==p2.and.q1/=q2) then
       ! we have two ray segments
       ya=(max(q1,q2)-1-vmod%cn)*vmod%dy+vmod%y0
       xa=linear_interp_1d(y1,x1,y2,x2,ya)
       ras(1,:)=(/x1,y1,xa,ya/)
       ras(2,:)=(/xa,ya,x2,y2/)
       pq(1,:)=(/p1,q1/)
       pq(2,:)=(/p2,q2/)
       nras=2
    end if

    if(p1/=p2.and.q1/=q2) then
       ! we have three ray segments
       ! we need to determine if we are first changiing cells in the x direction or y direction
       xa=(max(q1,q2)-1)*vmod%dy+vmod%y0-x1
       ya=(max(q1,q2)-1)*vmod%dy+vmod%y0-y1

       if (xa>=ya) then
          ! we are changing first in the y direction
          ya=(max(q1,q2)-1-vmod%cn)*vmod%dy+vmod%y0
          xa=linear_interp_1d(y1,x1,y2,x2,yb)
          xb=(max(p1,p2)-1-vmod%cn)*vmod%dx+vmod%x0
          yb=linear_interp_1d(x1,y1,x2,y2,xb)
          ras(1,:)=(/x1,y1,xa,ya/)
          ras(2,:)=(/xa,ya,xb,yb/)
          ras(3,:)=(/xb,yb,x2,y2/)
          pq(1,:)=(/p1,q1/)
          pq(2,:)=(/p1,q2/)
          pq(3,:)=(/p2,q2/)
          nras=3

       else
          ! we are changing first in the x direction
          xa=(max(p1,p2)-1-vmod%cn)*vmod%dx+vmod%x0
          ya=linear_interp_1d(x1,y1,x2,y2,xa)
          yb=(max(q1,q2)-1-vmod%cn)*vmod%dy+vmod%y0
          xb=linear_interp_1d(y1,x1,y2,x2,yb)
          ras(1,:)=(/x1,y1,xa,ya/)
          ras(2,:)=(/xa,ya,xb,yb/)
          ras(3,:)=(/xb,yb,x2,y2/)
          pq(1,:)=(/p1,q1/)
          pq(2,:)=(/p2,q1/)
          pq(3,:)=(/p2,q2/)
          nras=3
       end if

    end if




    !loop over all path segments

    do k=1,nras

       xa=ras(k,1)
       ya=ras(k,2)
       xb=ras(k,3)
       yb=ras(k,4)

       call cubic_bspline_basis(vmod,xa,ya,as,at)
       call cubic_bspline_interp_val(vmod,xa,ya,va)

       call cubic_bspline_basis(vmod,xb,yb,bs,bt)
       call cubic_bspline_interp_val(vmod,xb,yb,vb)


!      dl=   sqrt((xa-xb)**2+(ya-yb)**2)

! great_circle_distance
 dl=great_circle_distance(xa,xb,ya,yb)

!write(*,*) 'a',xa,ya
!write(*,*) 'b',xb,yb
!write(*,*) 'dl',dl
!if (isnan(dl)) then
!write(*,*) dl
!pause
!end if

       do i=1,4
          do j=1,4
             ii=pq(k,1)-2+i
             jj=pq(k,2)-2+j
             ra=as(i)*at(j)/va**2
             rb=bs(i)*bt(j)/vb**2
           
             if (ii<=vmod%nx+2*vmod%cn.and.jj<=vmod%ny+2*vmod%cn.and.ii>=1.and.jj>=1) then

  call velnod2id(ii,jj,vmod,nodid)

!write(*,*) dl
!write(*,*) ra,rb
                frero%val(nodid)=frero%val(nodid)-dl*(ra+rb)/2.0
                frero%nze(nodid)=1
             end if

          end do
       end do
    end do


  end subroutine calc_frechet_pathseg

end module m_arrivals



