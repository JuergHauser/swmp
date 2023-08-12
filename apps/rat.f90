! * rat *
!
! Program for multi arrival wavefront tracking and ray tracing for
! surface waves.
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

program rat

  use my_types
  use my_constants
  use my_functions

  use m_inout
  use m_solver
  use m_arrivals
  use m_wavefronts

  implicit none

  ! Data strcutures used
  type(scalar_field_2d)     :: vmod
  type(receivers)           :: recs
  type(receiver_mode)       :: recmode
  type(wavefront_container) :: wafc
  type(sources)             :: sous
  type(rat_conf)            :: conf
  type(node),pointer        :: head=>null()
  type(node),pointer        :: tail=>null()

  character(len=strlen)     :: arg

  ! counter variables
  integer:: h,i

  ! read in the configuration and allocate the relevant structures
  call get_command_argument(1,arg)
  call read_rat_conf(trim(arg),sous,recs,recmode,wafc,vmod,conf)

  ! setup frechet matrix files if needed
  if (conf%do_frechet==1) then
     call setup_frechet_matrix_files(vmod,conf)
  end if


  ! do the computations of the wavefronts for all source points and store
  ! the results in the wavefront container

  do h=1,sous%nn
     conf%sid=h ! sourcepoint index
     write(*,*) 'sourcepoint ',h
     write(*,*) 'wavefront tracking...'
     call reset_wavefronts_and_arrivals(wafc,recs,recmode,conf)
     call setup_strip_for_point_source(vmod,sous,recs,conf,head,tail)


     ! inital situation time t=0 gets stored as first iteration
     conf%wt%curit=1
     ! call receiver_record(head,recs,conf)
     call store_strip(head,tail,wafc,conf)

     ! start the computation i.e. evolution of the bicharateristic strip
     do i=2,conf%wt%maxit
        conf%wt%curit=i
       ! write(*,*) i,conf%wt%n_nod
        ! decide which ode solver to use
        select case (conf%wt%mode)

        case(1)
           call update_strip(head,vmod,conf)
        case(2)
           call update_strip(head,vmod,conf)
           call update_spfa(head,tail)
        end select


        ! A bicharacteristic strip with 3 points is describing wavefront,
        ! with only afew poitns in the computational domain.
        ! It therefore makes no sense to evolve it further
        ! and the computation for this source point is finished.

        call remove_nodes(head,tail,conf)

        if (conf%wt%n_nod<=3) then
           exit
        end if

        ! interpolating new nodes, if neccessary
        select case(conf%wt%interp)
        case(1)
           call insert_nodes_linear(head,tail,conf)
        case(2)
           call insert_nodes_splines(head,tail,conf)
        end select

        ! store the current computation step
        call store_strip(head,tail,wafc,conf)
        select case (conf%wt%mode)
        case(2)
           call update_dinr(head,tail)
        end select

     end do

     write(*,*) 'arrival information extraction...'
     call receiver_processing(recs,wafc,conf)
     if (conf%do_rays>0) then
        call extract_ray_paths_for_point_source(recs,wafc,conf)
        call analyze_ray_paths_for_point_source(recs,conf)
        ! compute the frechet matrix
        if (conf%do_frechet==1) then
           call save_frechet_matrix(recs,vmod,conf)
        end if
     end if

     ! start to write results to disk
     call write_wavefronts(wafc,conf)

     if (conf%out_rays>0) then
        call write_raypaths(recs,conf)
     end if

     call write_arrival_data(recs,conf%ofn_arrivals,conf%sid,conf%do_rays,conf%wt%mode)

  end do

  ! write summary
  call write_rat_summary(vmod,conf,recs)

  ! freeing up remaining memory
  call deallocate_wavefront_container(wafc)

  call deallocate_receivers(recs)
if (conf%recmode==1) then
  call deallocate_receiver_mode(recmode)
end if
  call deallocate_sources(sous)
  call delete_strip(head,tail)
  call deallocate_scalar_field_2d(vmod)
  call deallocate_rat_conf(conf)

end program rat
