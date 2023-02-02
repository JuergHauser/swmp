
module swmp


  use my_types
  use my_constants
  use my_functions

  use m_inout
  use m_solver
  use m_arrivals
  use m_wavefronts

  use iso_c_binding

  implicit none

! Data stuctures required for forward problem:
  type(scalar_field_2d)     :: vmod
  type(receivers)           :: recs
  type(receiver_mode)       :: recmode
  type(wavefront_container) :: wafc
  type(sources)             :: sous
  type(rat_conf)            :: conf
  type(node),pointer        :: head=>null()
  type(node),pointer        :: tail=>null()

! Data stuctures required for inverse problem:
  type(observed_arrivals):: obs
  type(predicted_arrivals):: pred

contains

! https://community.intel.com/t5/Intel-Fortran-Compiler/Passing-by-reference-a-Python-string-to-Fortran-function/m-p/1297543

subroutine read_conf(fn_ptr,fn_ptr_length)bind(c,name="rat_read_conf")
    type(c_ptr), value::  fn_ptr
    integer(c_int),value :: fn_ptr_length
    character(len=fn_ptr_length,kind=c_char), pointer :: fn_str

    ! local variables
    character(15)      :: ifn_vmod
    integer            :: i,j,k
    integer,pointer :: recsel(:)

    call c_f_pointer(fn_ptr, fn_str)
    call read_rat_conf(fn_str,sous,recs,recmode,wafc,vmod,conf)


end subroutine read_conf


subroutine forward()bind(c,name="forward")

 ! counter variables
  integer:: h,i

! setup frechet matrix files if needed
  if (conf%do_frechet==1) then
     call setup_frechet_matrix(vmod,conf)
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
           call calc_frechet_matrix(recs,vmod,conf)
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

end subroutine forward


! Setter and Getters to set everything from Python

! dt - ode solver time step size
  subroutine get_dt(dt)bind(c,name="get_dt")
  real(kind=c_float),intent(out) :: dt
    dt=conf%wt%dt
  end subroutine get_dt

  subroutine set_dt(dt)bind(c,name="set_dt")
    real(kind=c_float),intent(in),value :: dt
    conf%wt%dt=dt
  end subroutine set_dt

! maxit - maximum number of iterations
 subroutine get_maxit(maxit)bind(c,name="get_maxit")
  integer(kind=c_int),intent(out) :: maxit
    maxit=conf%wt%maxit
  end subroutine get_maxit

  subroutine set_maxit(maxit)bind(c,name="set_maxit")
    integer(kind=c_int),intent(in),value :: maxit
    conf%wt%maxit=maxit
  end subroutine set_maxit

  subroutine get_model_meta_data(x0,y0,nx,ny,dx,dy,cn) bind(c,name='get_model_meta_data')
    real(kind=c_float),intent(out) :: x0,y0,dx,dy
    integer(kind=c_int),intent(out) :: nx,ny,cn

    x0=vmod%x0
    y0=vmod%y0
    nx=vmod%nx
    ny=vmod%ny
    dx=vmod%dx
    dy=vmod%dy
    cn=vmod%cn

  end subroutine get_model_meta_data

  subroutine get_model_vector(val,nx,ny,cn)bind(c,name='get_model_vector')

    integer(kind=c_int), intent (in) :: nx,ny,cn
    real(c_float), intent(out) :: val((nx+cn*2)*(ny+cn*2))

    integer i,j,k

    do i=1,vmod%nx+vmod%cn*2
       do j=1,vmod%ny+vmod%cn*2
          call velnod2id(i,j,vmod,k)
          val(k)=vmod%val(i,j)
       end do
    end do

  end subroutine get_model_vector

  subroutine set_model_vector(val,nx,ny,cn)bind(c,name='set_model_vector')

    integer(kind=c_int), intent (in) :: nx,ny,cn
    real(c_float), intent(in) :: val((nx+cn*2)*(ny+cn*2))

    integer i,j,k

    do i=1,vmod%nx+vmod%cn*2
       do j=1,vmod%ny+vmod%cn*2
          call velnod2id(i,j,vmod,k)
          vmod%val(i,j)=val(k)
       end do
    end do

  end subroutine set_model_vector



  subroutine get_model_grid(val,nx,ny,cn)bind(c,name='get_model_grid')

    integer(kind=c_int), intent (in) :: nx,ny,cn
    real(c_float), intent(out) :: val((nx+cn*2),(ny+cn*2))

    integer i,j,k

    do i=1,vmod%nx+vmod%cn*2
       do j=1,vmod%ny+vmod%cn*2
          val(i,j)=vmod%val(i,j)
       end do
    end do

  end subroutine get_model_grid



  subroutine read_observations(fn_ptr,fn_ptr_length) bind(c,name='read_observations')

    type(c_ptr), value::  fn_ptr
    integer(c_int),value :: fn_ptr_length
    character(len=fn_ptr_length,kind=c_char), pointer :: fn_str
    call c_f_pointer(fn_ptr, fn_str)
    call read_observed_arrivals(obs,fn_str)
  end subroutine read_observations



  subroutine read_predictions(fn_ptr,fn_ptr_length) bind(c,name='read_predictions')
    type(c_ptr), value::  fn_ptr
    integer(c_int),value :: fn_ptr_length
    character(len=fn_ptr_length,kind=c_char), pointer :: fn_str
    call c_f_pointer(fn_ptr, fn_str)
    call read_predicted_arrivals(pred,fn_str)
  end subroutine read_predictions


 subroutine get_number_of_observations(n) bind(c,name='get_number_of_observations')
     integer(kind=c_int), intent(out) :: n
     n =obs%n
 end subroutine get_number_of_observations

  subroutine get_obervations(val,n) bind(c,name='get_observations')
    integer(kind=c_int), intent(out) :: n
    real(c_float), intent(out) :: val(obs%n,6)
    integer i

  do i=1,obs%n
    val(i,1) =obs%sou(i)
    val(i,2) =obs%rec(i)
    val(i,3) =obs%arn(i)
    val(i,4) =obs%art(i)
    val(i,5) =obs%azi(i)
    val(i,6) =obs%unc(i)
 end do


end subroutine get_obervations


    subroutine get_number_of_predictions(n) bind(c,name='get_number_of_predictions')
     integer(kind=c_int), intent(out) :: n
     n =pred%n
    end subroutine get_number_of_predictions


  subroutine get_predictions(val,n) bind(c,name='get_predictions')
  integer(kind=c_int), intent(in) :: n
    real(c_float), intent(out) :: val(pred%n,5)
    integer i

    do i=1,pred%n
       val(i,1) = pred%sou(i)
       val(i,2) = pred%rec(i)
       val(i,3) = pred%arn(i)
       val(i,4) = pred%art(i)
       val(i,5) = pred%azi(i)
    end do

  end subroutine get_predictions




end module swmp


module modgen

  use my_functions
  use my_types
  use m_inout

    use iso_c_binding

  implicit none
  type(scalar_field_2d) :: vel
  type(scalar_field_2d) :: cov
  type(scalar_field_2d)	:: ani
  type(gen2dv_conf)     :: conf
  type(list_ind_val)    :: spikes

contains

subroutine read_conf(fn_ptr,fn_ptr_length)bind(c,name="gen2dv_read_conf")
    type(c_ptr), value::  fn_ptr
    integer(c_int),value :: fn_ptr_length
    character(len=fn_ptr_length,kind=c_char), pointer :: fn_str

    ! local variables
    character(15)      :: ifn_vmod
    integer            :: i,j,k
    integer,pointer :: recsel(:)

    call c_f_pointer(fn_ptr, fn_str)
    call read_gen2dv_conf(fn_str,conf,spikes)


end subroutine read_conf

subroutine generate() bind(c,name="gen2dv_run")


  real(kind=dbl)::chvp,chvp1,chvp2
  integer ::vusp1,vusp2,vusp1o,vusp2o
  integer i,j


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



end subroutine generate


end module modgen



