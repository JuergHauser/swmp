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
! * my_types *
! * my_constants *
!
! user derived types and constants used throughout swmp
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

module my_types
  ! this module provides the user derived types
  use, intrinsic :: ISO_C_BINDING
  implicit none

  !****************************************************************************!
  ! parameter declarations                                                     !
  !****************************************************************************!

  integer, parameter :: sgl = selected_real_kind(p=6,r=37)
  integer, parameter :: dbl = selected_real_kind(p=13,r=200)
  integer, parameter :: strlen = 64

  integer,parameter :: input=101
  integer,parameter :: output=102

  integer,parameter :: bin=103
  integer,parameter :: hdr=104

  !****************************************************************************!
  ! dummy variables                                                            !
  !****************************************************************************!

  character(strlen) :: cdum
  integer           :: idum
  real(kind=dbl)    :: rdum


  !****************************************************************************!
  ! types for vectors, matrices and lists                                      !
  !****************************************************************************!

  type scalar_field_2d
     real(kind=dbl),pointer :: val(:,:)=>null() ! values
     real(kind=dbl)         :: x0,y0            ! lower left corner
     real(kind=dbl)         :: dx,dy            ! grid spacing
     integer                :: nx,ny            ! number of grid points nx*ny
     integer                :: cn               ! number of cushion nodes
  end type scalar_field_2d

  type matrix
     real(kind=dbl),pointer  :: val(:,:)=>null() ! values
     integer                 :: n1,n2            ! number of elements n1xn2
  end type matrix

  type vector
     real(kind=dbl),pointer :: val(:)=> null()  !
     integer :: n                               !
  end type vector

 type list
     integer,pointer :: val(:)=> null()  !
     integer :: n
  end type list

  type crs_matrix                               ! compact row storage matrix
     integer,pointer         :: col(:)=>null()
     integer,pointer         :: row(:)=>null()
     real(kind=dbl),pointer  :: val(:)=>null()
     integer :: n1,n2
  end type crs_matrix

  type list_index
     integer            :: n
     integer,pointer    :: ind(:,:) => null()
  end type list_index




  type list_ind_val
     integer                 :: n                ! number of nodes with
     integer,pointer         :: ind(:,:)=>null() ! index
     real(kind=dbl),pointer  :: val(:)=>null()   ! value
  end type list_ind_val



  !****************************************************************************!
  ! configuration files                                                        !
  !****************************************************************************!

  type wavefront_tracker_configuration
     real(kind=dbl)  :: dt
     integer         :: solv
     integer         :: mode
     integer         :: maxit
     integer         :: curit
     integer         :: ipmod
     integer         :: n_senod
     integer         :: n_nod
     real(kind=dbl)  :: pdmi,pdma
     real(kind=dbl)  :: scax,scay
     integer         :: wafso(2)
     integer         :: nfrid
     integer         :: frech_do
     integer         :: rapex_do
     integer         :: init_nodes
     integer          :: solver
     integer          :: interp
   !  complex(kind=dbl):: epsilon
  end type wavefront_tracker_configuration

  type creobs_conf
     character(len=strlen) :: ifn_traveltimes
     character(len=strlen) :: ofn_traveltimes
     integer :: do_ran
     real(kind=dbl):: stdev
     integer :: rseed
     real(kind=dbl):: unc
  end type creobs_conf

  type gcray_conf
     character(len=strlen)      :: ifn_sources
     character(len=strlen)      :: ifn_receivers
     character(len=strlen)      :: ofn_traveltimes
     real(kind=dbl):: vel
          real(kind=dbl):: rearth

  end type gcray_conf

  type gen2dv_conf  ! params for gen2dv.f90
     real(kind=dbl) :: x0,y0     ! lower left corner
     integer        :: nx,ny     ! number of grid cells
     real(kind=dbl) :: dx,dy     ! grid spacing
     integer        :: cn        ! number of cushion nodes
     real(kind=dbl) :: vbg       ! background velocity
     integer        :: do_modcov  ! add model covariance
     real(kind=dbl) :: covbg       ! background covariance
     integer        :: do_rand  ! add random structure
     real(kind=dbl) :: stadev    ! standard deviation of random stucture
     integer        :: rseed     ! random seed for Gaussian noise
     integer        :: do_cebo    ! do checker board
     real(kind=dbl)       :: cebo_pert ! do specific spikes
     integer      :: cebo_size
     integer   :: cebo_spac
     character(len=strlen):: ofn_velmod
     character(len=strlen):: ofn_modcov
     character(len=strlen):: ofn_animod
  end type gen2dv_conf



  type rat_conf
     character(len=strlen)      :: ifn_velmod
     character(len=strlen)      :: ifn_sources
     character(len=strlen)      :: ifn_receivers
     character(len=strlen):: ifn_recmode
     character(len=strlen)      :: ofn_arrivals
     character(len=strlen)      :: ofn_wavefronts
     character(len=strlen)      :: ofn_raypaths
     character(len=strlen)      :: ofn_frechet_hdr
     character(len=strlen)      :: ofn_frechet_mat
     character(len=strlen)      :: ofn_frechet_rai
     character(len=strlen)      :: ofn_sumfile
     type(wavefront_tracker_configuration) :: wt
     real(kind=dbl):: xmin,xmax,ymin,ymax
     integer :: out_wafint
     integer :: do_rays
     integer :: out_rays
     integer :: do_frechet
     integer,pointer    :: tonar(:)=>null()
     integer :: velint
     integer :: sid
          real(kind=dbl):: rearth

integer :: recmode
     integer :: enable_file_output  ! 0=no disk writes, 1=write files (for backward compat)
  end type rat_conf

  type inv_conf
     integer :: subdim					! subspace dimensions
     integer :: invmod					! inversion mode
     real(kind=dbl):: svdthre ! svd threshold for singualr values..
     character(len=strlen)     :: ifn_vinit          ! inital velocity model
     character(len=strlen)     :: ifn_vcur           ! current velocity model
     character(len=strlen)	  :: ifn_vunc
     character(len=strlen)     :: ifn_predtimes        ! predicted arrrival times
     character(len=strlen)     :: ifn_obstimes       ! observed arrival times
     character(len=strlen)     :: ifn_frechet_hdr        ! freche matrix
     character(len=strlen)     :: ifn_frechet_mat        ! freche matrix
     character(len=strlen)     :: ifn_frechet_rai    !
     character(len=strlen)     :: ifn_info           ! rat.inof file
     character(len=strlen)     :: ofn_vpert          ! velcoity perturbation
     character(len=strlen)     :: ofn_vupd           ! update velocity field
     character(len=strlen)     :: ofn_tres           ! traveltime_residuals
     character(len=strlen)     :: ofn_rapco          ! raypath coverage map
     integer                   :: do_sds              ! secodn derivative smoothing
     real(kind=dbl):: rearth

  end type inv_conf


type rabgv_conf
     character(len=strlen):: ifn_mod
     character(len=strlen):: ifn_rac
     real(kind=dbl):: bgv
real(kind=dbl):: mnr
 character(len=strlen):: ofn_mod
 
end type rabgv_conf

  type sia_conf
     character(len=strlen):: ifn_vmoda
     character(len=strlen):: ifn_vmodb
     character(len=strlen):: ifn_atimes
     character(len=strlen):: ifn_btimes
     character(len=strlen):: ifn_obstimes
     integer :: lk,ak,u
     integer :: lmax,amin,umax
     real(kind=dbl):: dv
     real(kind=dbl):: omegamin,omegamax
     integer :: smof_type,smof_size
     character(len=strlen):: ofn_vmoda
     character(len=strlen):: ofn_atimes
     character(len=strlen):: ofn_vmodb
     character(len=strlen):: ofn_sumfile
     integer :: spert,sacep
     real(kind=dbl):: tk
     integer :: k
  end type sia_conf

  type per2dv_conf
     character(len=strlen):: ifn_vmod
     character(len=strlen):: ofn_vmod
     real(kind=dbl):: dv
     integer :: seed
  end type per2dv_conf

  type tresid_conf
     character(len=strlen)::ifn_predtimes
     character(len=strlen)::ifn_obstimes
     integer              :: dump_res
     character(len=strlen):: ofn_residuals
     character(len=strlen):: ofn_sumfile
  end type tresid_conf


  type resa2d_conf    ! params for resa2dv.f90
     character(len=strlen) :: igrid            ! filename
     character(len=strlen) :: ovef,ogaf        ! filename
     character(len=strlen) :: gmtvef,gmtgaf    ! filenmae
     integer       :: nrx, nry         ! resmapling factor
     integer       :: sw               ! switch
  end type resa2d_conf


type modro_conf
     character(len=strlen) :: ifn_vcur
     character(len=strlen) :: ifn_vref
     integer               :: nrx,nry
     real(kind=dbl):: rearth

end type modro_conf

  type misfit_conf
     character(len=strlen) :: ifn_vref     ! referenz model
     character(len=strlen) :: ifn_vinv     ! inverted model
     character(len=strlen) :: ifn_vtrue     ! true model
     character(len=strlen) :: ofn_dabs     ! absoulte diefference
     character(len=strlen) :: ofn_drel    ! relative difference
     character(len=strlen) :: ofn_adabs     ! absoulte anomaly difference
     character(len=strlen) :: ofn_adrel    ! relative anomaly difference
     character(len=strlen):: ofn_sumfile
     integer               :: nrx, nry         ! resmapling factor
  end type misfit_conf

type vn2vj_conf
     character(len=strlen) :: ifn_vmod
     integer :: cn
     character(len=strlen) :: ofn_vmod
end type vn2vj_conf

  !****************************************************************************!
  ! specific files for wavefront tracking                                      !
  !****************************************************************************!

  type raypath
     integer :: n
     real(kind=dbl),pointer :: pos(:,:)=>null()
  end type raypath


  type arrival
     real(kind=dbl)    :: time   ! arrival time
     integer           :: iter   ! iteration
     real(kind=dbl)    :: azi    ! angle of incidence
     real(kind=dbl) :: spf       ! spreading factor
     integer(kind=dbl) :: raid   ! ray id
     type(raypath)     :: path   ! raypath
     integer		   :: paic   ! path in coputationa domain
     real(kind=dbl):: gaba ! gaussian bema amplitude
  end type arrival

  type sources
     integer :: type     ! (1) point source (2) plane wave
     real (kind=dbl),pointer :: pos(:,:)=>null() ! source position
     real(kind=dbl),pointer  :: agi(:)=>null()   ! angle
     integer                 :: nn               ! number of source points
  end type sources

  type receivers
     real(kind=dbl),pointer :: pos(:,:)=>null()       ! position of the receiver
     integer                :: nn                     ! number of receviers
     integer                :: mar                    !  maximum number of arrivals
     integer,pointer        :: nar(:)=>null()         ! number of arrivals
     type(arrival),pointer	:: ariv(:,:)
integer,pointer :: stat(:)
  end type receivers

  type receiver_mode
     integer :: nn
     type(list),pointer:: act(:)
  end type receiver_mode

  type wavefront
     integer                :: headpid
     integer                :: n_nod
     real(kind=dbl),pointer :: pos(:,:)

     real(kind=dbl),pointer :: spf(:)
     integer,pointer        :: neig(:,:)
     integer,pointer        :: con(:,:)
     integer,pointer        :: stat(:)
     type(list_index)       :: index
  end type wavefront

  type wavefront_container
     type(wavefront),pointer ::waf(:)
     integer :: n
  end type wavefront_container


  type node
     type(node),pointer  :: prev=>null()   ! pointer to the prev neighbours
     type(node),pointer  :: next=>null()   ! pointer to then next neighbours
     real(kind=dbl)      :: pos(3)         ! position
     real(kind=dbl)      :: spf
     real(kind=dbl)      :: dinr
     integer             :: con(2)         ! connected
     integer             :: pid            ! unique point id
  end type node


  !****************************************************************************!
  ! specific files for the inversion                                           !
  !****************************************************************************!

  type observed_arrivals
     integer :: n
     real(kind=dbl),pointer::art(:)
     integer,pointer :: sou(:)
     integer,pointer :: rec(:)
     integer,pointer :: arn(:)
     real(kind=dbl),pointer :: unc(:)
     real(kind=dbl),pointer :: azi(:)
  end type observed_arrivals

  type predicted_arrivals
     integer :: n
     real(kind=dbl),pointer::art(:)
     integer,pointer :: sou(:)
     integer,pointer :: rec(:)
     integer,pointer :: arn(:)
     real(kind=dbl),pointer :: azi(:)
  end type predicted_arrivals

  type traveltime_residuals
     integer :: n
     integer :: tnar
     real(kind=dbl):: trms
     real(kind=dbl):: tvar
     real(kind=dbl):: tmean
     type(vector) :: rms
     type(vector) :: var
     type(vector) :: mean
     type(vector) :: nar
     type(vector),allocatable:: ares(:)
  end type traveltime_residuals


  type frechet_matrix_header
     integer :: nre                 ! number of records
     integer :: n1,n2
     integer :: recole
  end type frechet_matrix_header

  type frechet_matrix_row
     integer,pointer :: nze(:)
     real(kind=dbl),pointer  :: val(:)
  end type frechet_matrix_row

contains

  ! these subroutine allow the deallocation of the pointers in the user derived types

  !----------------------------------------------------------------------------!

  subroutine  deallocate_scalar_field_2d(a)
    implicit none
    type(scalar_field_2d)::a
    deallocate(a%val);nullify(a%val)
    a%x0=0.0;a%y0=0.0
    a%nx=0;a%ny=0;a%cn=0
    a%dx=0;a%dy=0
  end subroutine deallocate_scalar_field_2d

  subroutine deallocate_matrix(a)
    implicit none
    type(matrix) :: a
    deallocate(a%val);nullify(a%val)
    a%n1=0;a%n2=0
  end subroutine deallocate_matrix

  subroutine deallocate_vector(a)
    implicit none
    type(vector)::a
    deallocate(a%val);nullify(a%val)
    a%n=0
  end subroutine deallocate_vector


subroutine deallocate_list(a)
implicit none
type(list):: a
 deallocate(a%val);nullify(a%val)
a%n=0
end subroutine deallocate_list


  subroutine deallocate_crs_matrix(a)
    implicit none
    type(crs_matrix):: a
    deallocate(a%val);nullify(a%val)
    deallocate(a%col);nullify(a%col)
    deallocate(a%row);nullify(a%row)
    a%n1=0;a%n2=0
  end subroutine deallocate_crs_matrix

  subroutine deallocate_sources(sous)
    implicit none
    type(sources) :: sous

    if ( associated(sous%pos)) then
       deallocate(sous%pos);nullify(sous%pos)
    end if

    if (associated(sous%agi)) then
       deallocate(sous%agi);nullify(sous%agi)
    end if

    sous%nn=0
  end subroutine deallocate_sources

  subroutine deallocate_receivers(recs)
    implicit none
    type(receivers) :: recs
    integer :: i,j
    do i=1,recs%nn
       do j=1,recs%nar(i)
          if (recs%ariv(i,j)%path%n/=0) then
             recs%ariv(i,j)%path%n=0
             deallocate (recs%ariv(i,j)%path%pos)
          end if
       end do
    end do
    deallocate(recs%ariv)
    deallocate(recs%pos)
    deallocate(recs%nar)
  end subroutine deallocate_receivers


subroutine deallocate_receiver_mode(recmode)
implicit none
type(receiver_mode):: recmode
integer ::i

do i=1,recmode%nn
call deallocate_list(recmode%act(i))
end do
deallocate(recmode%act)
recmode%nn=0
end subroutine deallocate_receiver_mode

  subroutine deallocate_rat_conf(conf)
    implicit none
    type(rat_conf):: conf
    deallocate(conf%tonar)
  end subroutine deallocate_rat_conf

  subroutine deallocate_wavefront_container(wafc)
    implicit none
    type(wavefront_container) ::wafc
    integer :: i
    do i=1,wafc%n
       if (wafc%waf(i)%n_nod>0) then

          deallocate(wafc%waf(i)%pos)
          deallocate(wafc%waf(i)%neig)
          deallocate(wafc%waf(i)%con)
          deallocate(wafc%waf(i)%stat)
          deallocate (wafc%waf(i)%index%ind)
          wafc%waf(i)%index%n=0
          wafc%waf(i)%n_nod=0

       end if
    end do
    deallocate(wafc%waf)
  end subroutine deallocate_wavefront_container

  subroutine deallocate_traveltime_residuals(tresid)
    type(traveltime_residuals):: tresid
    integer :: i
    call deallocate_vector(tresid%rms)
    call deallocate_vector(tresid%nar)
    do i=1,tresid%n
       call deallocate_vector(tresid%ares(i))
    end do
    deallocate(tresid%ares)

  end subroutine deallocate_traveltime_residuals


 subroutine deallocate_predicted_arrivals(pred)
type(predicted_arrivals):: pred
pred%n=0
deallocate(pred%art)
deallocate(pred%sou)
deallocate(pred%rec)
deallocate(pred%arn)
deallocate(pred%azi)
   end subroutine deallocate_predicted_arrivals


 subroutine deallocate_observed_arrivals(obs)
type(observed_arrivals):: obs
obs%n=0
deallocate(obs%art)
deallocate(obs%sou)
deallocate(obs%rec)
deallocate(obs%arn)
deallocate(obs%azi)
deallocate(obs%unc)
   end subroutine deallocate_observed_arrivals


  subroutine deallocate_list_index(list)
    implicit none
    type(list_index) :: list
    deallocate(list%ind);nullify(list%ind)
  end subroutine deallocate_list_index

  !------------------------------------------------------------------------------!

  subroutine reallocate_list_index(r1)
    ! reallocates a list(roster) of indices
    implicit none
    ! subroutien arguments
    type(list_index):: r1
    ! local variables
    type(list_index):: r2
    integer ::n


    n=size(r1%ind,2)

    r2%n=r1%n
    allocate(r2%ind(r2%n,n))
    r2%ind=r1%ind(1:r1%n,:)
    call deallocate_list_index(r1)
    r1%n=r2%n
    allocate(r1%ind(r1%n,n))
    r1%ind=r2%ind

    call deallocate_list_index(r2)

  end subroutine reallocate_list_index


  subroutine reallocate_vector(v1)
    ! reallocates a vector
    implicit none
    ! subroutien arguments
    type(vector):: v1
    ! local variables
    type(vector):: v2

    v2%n=v1%n
    allocate(v2%val(v2%n))
    v2%val=v1%val(1:v1%n)
    call deallocate_vector(v1)
    v1%n=v2%n
    allocate(v1%val(v1%n))
    v1%val=v2%val

  end subroutine reallocate_vector

  !------------------------------------------------------------------------------!

end module my_types

!=============================================================================!

module my_constants
  ! This module provides the constants pi and e as given by Mathematica
  use my_types
  implicit none
  real(kind=dbl) :: rearth=6371.0_dbl
  real(kind=dbl),parameter :: radian=180.0_dbl/3.14159265358979323846264338327950288419716939937
  real(kind=dbl),parameter :: pi=3.14159265358979323846264338327950288419716939937
  real(kind=dbl),parameter :: euler =2.71828182845904523536028747135266249775724709369
end module my_constants
