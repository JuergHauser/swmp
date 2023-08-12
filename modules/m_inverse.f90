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
! * m_inverse *
!
! svd subspace etc
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

module m_inverse

  implicit none
contains



!-----------------------------------------------------------------------------!


  subroutine compute_residuals(obs,pred,resid)

    use my_types

    implicit none

    ! subroutine arguments
    type(observed_arrivals)  :: obs
    type(predicted_arrivals) :: pred
    type(traveltime_residuals):: resid

    ! local variables
    type(list_index)          :: usarr
    integer :: i,j
    integer :: arn


    resid%tnar=0
    resid%trms=0
    resid%tmean=0
   resid%tvar=0


    usarr%n=maxval((/obs%n,pred%n/))
    allocate(usarr%ind(usarr%n,2))
    usarr%n=0

    do i=1,obs%n
       do j=1,pred%n
          if (obs%sou(i)==pred%sou(j).and.obs%rec(i)==pred%rec(j).and.obs%arn(i)==pred%arn(j)) then
             usarr%n=usarr%n+1
             usarr%ind(usarr%n,:)=(/i,j/)
          end if
       end do
    end do


    ! The arrivals are matchhed so now it has to be determined how many arrivals will be used
    ! That means classes of later arrivals.

    resid%n=max(maxval(obs%arn(:)),maxval(pred%arn(:)))

    ! allocating the storage

    resid%rms%n=resid%n
    allocate(resid%rms%val(resid%rms%n))
   resid%mean%n=resid%n
    allocate(resid%mean%val(resid%rms%n))
   resid%var%n=resid%n
    allocate(resid%var%val(resid%rms%n))

    resid%nar%n=resid%n
    allocate(resid%nar%val(resid%nar%n))
    resid%nar%val=0.0_dbl
    allocate(resid%ares(resid%n))

    ! count the number of arrivals

    do i=1,usarr%n
       arn=obs%arn(usarr%ind(i,1))
       resid%nar%val(arn)=resid%nar%val(arn)+1
    end do

    ! allocate the storage for the residuals
    do i=1,resid%n
       resid%ares(i)%n=resid%nar%val(i)
       allocate(resid%ares(i)%val(resid%ares(i)%n))
       resid%ares(i)%n=0
    end do



    !  working on the different arrival classes...
    do i=1,usarr%n
       arn=obs%arn(usarr%ind(i,1))
       resid%ares(arn)%n=resid%ares(arn)%n+1
       resid%ares(arn)%val(resid%ares(arn)%n)=pred%art(usarr%ind(i,2))-obs%art(usarr%ind(i,1))
    end do



    do i=1,resid%n
       if (resid%nar%val(i)>0) then
          resid%rms%val(i)=sqrt(sum((resid%ares(i)%val)**2)/resid%nar%val(i))
          resid%mean%val(i)=sum(resid%ares(i)%val)/resid%nar%val(i)
          resid%var%val(i)=sum(resid%ares(i)%val-resid%mean%val(i))**2/resid%nar%val(i)
       else
          resid%rms%val(i)=0.0
          resid%mean%val(i)=0.0
          resid%var%val(i)=0.0
       end if
    end do


    ! compute trms values
    do i=1,resid%n
       resid%tnar=resid%tnar+resid%nar%val(i)
       resid%trms=resid%trms+sum((resid%ares(i)%val)**2)
       resid%tmean=resid%tmean+sum(resid%ares(i)%val)
    end do

    resid%trms=sqrt(resid%trms/resid%tnar)
    resid%tmean=resid%tmean/resid%tnar

    do i=1,resid%n
resid%tvar=resid%tvar+sum((resid%ares(i)%val-resid%tmean)**2)
    end do
resid%tvar=resid%tvar/resid%tnar

  end subroutine compute_residuals

!-----------------------------------------------------------------------------!

subroutine initialise_data(obs,pred,dobs,cdi,dpred,usarr)

  use my_types
  implicit none


  ! subroutine arguments
  type(observed_arrivals)   :: obs
  type(predicted_arrivals)  :: pred
  type(vector)              :: dobs    ! observation
  type(vector)              :: dpred   ! model prediction
  type(list_index)          :: usarr
  type(crs_matrix):: cdi

  ! local variables
  integer :: i,j,k


  usarr%n=maxval((/obs%n,pred%n/))
  allocate(usarr%ind(usarr%n,2))
  usarr%n=0


  do i=1,obs%n
     do j=1,pred%n
        if (obs%sou(i)==pred%sou(j).and.obs%rec(i)==pred%rec(j).and.obs%arn(i)==pred%arn(j)) then
           usarr%n=usarr%n+1
           usarr%ind(usarr%n,:)=(/i,j/)
        end if
     end do
  end do

  ! usarr is now the list of arrivals which will be used...

  call reallocate_list_index(usarr)

  ! observations with error

  dobs%n=usarr%n
  allocate(dobs%val(dobs%n))
  do i=1,usarr%n
     dobs%val(i)=obs%art(usarr%ind(i,1))
  end do


  ! The data covariance matrix cd is adiagonal matrix
  ! As we only use (cd)^-1 we can store the ivnerse directly i.e. cdi


  cdi%n1=usarr%n
  cdi%n2=usarr%n

  allocate(cdi%val(cdi%n1))
  allocate(cdi%col(cdi%n2))
  allocate(cdi%row(cdi%n1+1))

  do i=1,usarr%n
     cdi%val(i)=1.0/obs%unc(usarr%ind(i,1))
     cdi%col(i)=i
     cdi%row(i)=i
  end do
  cdi%row(cdi%n1+1)=usarr%n+1

  ! predictions

  dpred%n=usarr%n
  allocate(dpred%val(dpred%n))
  do j=1,usarr%n
     dpred%val(j)=pred%art(usarr%ind(j,2))
  end do


end subroutine initialise_data


  !----------------------------------------------------------------------------!

  subroutine initialise_model(m0,vcur,m,cm,cmi,conf)

    use my_types
    use my_functions
    use m_inout
use m_linalg
    implicit none

    ! subroutine arguments
    type(scalar_field_2d) :: vcur
    type(crs_matrix):: cm,cmi
    type(inv_conf):: conf
    type(vector) :: m, m0

    ! local variables
    integer :: i,j,k
    type(scalar_field_2d):: vunc,vinit


    ! load vinit and rewrite into m0
    call read_scalar_field_2d(vinit,conf%ifn_vinit)
    m0%n=(vinit%nx+vinit%cn*2)*(vinit%ny+vinit%cn*2)
    allocate(m0%val(m0%n))
    do i=1,vinit%nx+vinit%cn*2
       do j=1,vinit%ny+vinit%cn*2
          call velnod2id(i,j,vinit,k)
          m0%val(k)=vinit%val(i,j)
       end do
    end do

    ! we no longer need vinit
    call deallocate_scalar_field_2d(vinit)

    ! load  vcur and rewrite into m
    call read_scalar_field_2d(vcur,conf%ifn_vcur)
    m%n=(vcur%nx+vcur%cn*2)*(vcur%ny+vcur%cn*2)
    allocate(m%val(m%n))
    do i=1,vcur%nx+vcur%cn*2
       do j=1,vcur%ny+vcur%cn*2
          call velnod2id(i,j,vcur,k)
          m%val(k)=vcur%val(i,j)
       end do
    end do

! load vunc and rewrite into cm
    call read_scalar_field_2d(vunc,conf%ifn_vunc)

    ! the model covariance matrix is a square diagonal matrix

    cm%n1=(vunc%nx+vunc%cn*2)*(vunc%ny+vunc%cn*2)
    cm%n2=cm%n1

    allocate(cm%col(cm%n1))
    allocate(cm%row(cm%n2+1))
    allocate(cm%val(cm%n1))

    k=0
    do i=1,vunc%nx+vunc%cn*2
       do j=1,vunc%ny+vunc%cn*2
          k=k+1
          cm%col(k)=k
          cm%row(k)=k
          cm%val(k)=vunc%val(i,j)
       end do
    end do

    cm%row(k+1)=cm%n1+1


    cmi=crsinversdiag(cm)



    ! we no longer need vunc
    call deallocate_scalar_field_2d(vunc)

  end subroutine initialise_model


  !----------------------------------------------------------------------------!

  subroutine compute_smoothing_operator(vcur,m,dtd,conf)

    use my_types
    use my_functions
    use m_linalg
    implicit none

    ! subroutine arguments
    type(crs_matrix) :: dtd
    type(scalar_field_2d):: vcur
    type(vector):: m
    type(inv_conf) :: conf

    ! local variables
    type(list_index):: nei
    type(crs_matrix) :: d,dt
    integer :: h,i,j,k
    integer :: nze

    type(matrix):: tmp


    ! initialise the neighbour list
    nei%n=0
    allocate(nei%ind(5,1)) ! there can't be more than four neighbours...

    ! allocate d

    ! we have to count the number of non zero elemnts in the
    ! smoothing operator D...

    nze=((vcur%nx+vcur%cn*2)-2)*((vcur%ny+vcur%cn*2)-2)*5 &
         +(vcur%nx+vcur%cn*2-2)*2*4+(vcur%ny+vcur%cn*2-2)*2*4 &
         +4*3

    ! we build D row by row

    d%n1=m%n
    d%n2=m%n

    allocate(d%val(nze))
    allocate(d%col(nze))
    allocate(d%row(m%n+1))
    d%val=0.0





    h=0
    k=1

    do i=1,m%n

       ! working for m(i)


       call id2nei(i,vcur,nei)
       nei%n=nei%n+1
       nei%ind(nei%n,1)=i
       call heap_sort_list_index(nei)

       do j=1,nei%n

          if (nei%ind(j,1)==i) then
             d%val(h+j)=-(nei%n-1)
             d%col(h+j)=nei%ind(j,1)
          else
             d%val(h+j)=1
             d%col(h+j)=nei%ind(j,1)
          end if
       end do
       h=h+nei%n
       d%row(i)=k
       k=k+nei%n
    end do

    d%row(m%n+1)=h+1

dt=crstranspose(d)
dtd=crsmatmul(dt,d)


call deallocate_crs_matrix(d)
call deallocate_crs_matrix(dt)


! we always build the smoothing oeprator dtd
! It is set to zero if conf$do_sds==0

if (conf%do_sds==0) then
dtd%val=0
end if

call deallocate_list_index(nei)

  end subroutine compute_smoothing_operator

!--------------------------------------------------------------------------------!

  subroutine compute_projection_matrix(g,gt,dobs,cdi,dpred,m0,cm,cmi,m,dtd,epsilon,eta,gamma,a,conf)
    ! This subroutien creates the projection matrix a using the hessian h
    use my_types
    use my_functions
    use m_linalg

    !subroutine arguments
    type(crs_matrix)::g,gt,cdi,cm,cmi,dtd
    type(vector):: dpred,dobs,m0,m,gamma
    real(kind=dbl):: epsilon,eta
    type(inv_conf):: conf
    type(matrix):: a


    ! local variables
    integer :: i,j,k
    real(kind=dbl):: thre
    type(crs_matrix):: mh
    type(vector)::mgamma

    ! some temporary storage gets only allocated when needed
    type(crs_matrix):: h
    type(crs_matrix) :: ctmp1,ctmp2
    type(crs_matrix) :: ctmpa,ctmpb,ctmpc
    type(matrix) :: mtmp1
    type(vector) :: vtmp1,vtmp2
    type(vector) :: vtmpa,vtmpb,vtmpc

    thre=0.0_dbl

    ! first we need to build gamma

    ! gamma= G^T*Cd^-1 (dpre-dobs) +epsilond Cm^-1*(m-m0) +eta*D^T*D*m

    ! G^T Cd^-1 (dpre-dobs)
    vtmp1=vecsub(dpred,dobs)
    gt=crstranspose(g)
    ctmp1=crsmatmul(gt,cdi)
    vtmpa=crsmatvecmul(ctmp1,vtmp1)


    call deallocate_vector(vtmp1)
    call deallocate_crs_matrix(ctmp1)

    ! epsilon Cm^-1*(m-m0)
    vtmp1=vecsub(m,m0)
    vtmp2=crsmatvecmul(cmi,vtmp1)
    vtmpb=scavecmul(epsilon,vtmp2)
    call deallocate_vector(vtmp1)
    call deallocate_vector(vtmp2)

    ! eta*D^TD*m
    ctmp1=crsscamatmul(eta,dtd)
    vtmpc=crsmatvecmul(ctmp1,m)
    call deallocate_crs_matrix(ctmp1)



    ! build gamma
    vtmp1=vecadd(vtmpa,vtmpb)
    gamma=vecadd(vtmp1,vtmpc)
    call deallocate_vector(vtmp1)
    call deallocate_vector(vtmpa)
    call deallocate_vector(vtmpb)
    call deallocate_vector(vtmpc)




    ! we need to build the hessian h
    ! h=G^T*Cd^-1*G + epsilon* Cm^-1 + eta*D^T*D

    ! G^T*Cd^-1*G

    ctmp1=crsmatmul(gt,cdi)
    ctmpa=crsmatmul(ctmp1,g)
    call deallocate_crs_matrix(ctmp1)

    ! epsilon* Cm^-1
    ctmpb=crsscamatmul(epsilon,cmi)

    ! eta*D^T*D
    ctmpc=crsscamatmul(eta,dtd)

    ! build hessian
    ctmp1=crsmatadd(ctmpa,ctmpb)
    h=crsmatadd(ctmp1,ctmpc)

    call deallocate_crs_matrix(ctmp1)
    call deallocate_crs_matrix(ctmpa)
    call deallocate_crs_matrix(ctmpb)
    call deallocate_crs_matrix(ctmpc)



    ! start to build the projection matrix a
    ! we can not store matrix a in crs storage. There is a risk
    ! that a contains a column of zeros these denies the possibility
    ! to store a as a crs matrix


    mgamma=crsmatvecmul(cm,gamma)

    mtmp1%n1=mgamma%n  ! number of model parameters
    mtmp1%n2=conf%subdim ! number of subdimenisons
    allocate(mtmp1%val(mtmp1%n1,mtmp1%n2))

    call normalize_vector(mgamma%val)

    vtmp1%n=mgamma%n
    allocate(vtmp1%val(vtmp1%n))
    vtmp1%val=mgamma%val

call deallocate_vector(mgamma)

    mtmp1%val(:,1)=vtmp1%val


    mh=crsmatmul(cm,h)

    do i=2,conf%subdim,1

       vtmp2=crsmatvecmul(mh,vtmp1)
 call normalize_vector(vtmp2%val)
       ! write(*,*) vtmp2%val

       mtmp1%val(:,i)=vtmp2%val

       vtmp1%n=vtmp2%n
       vtmp1%val=vtmp2%val
       call deallocate_vector(vtmp2)

    end do


    call deallocate_vector(vtmp1)

    !othogonalize a
    a=svd_orthogonalization(mtmp1,conf%svdthre)

    call deallocate_matrix(mtmp1)

    ! deallocate the hessian matrix, as it is no longer needed
    call deallocate_crs_matrix(h)
    call deallocate_crs_matrix(mh)

  end subroutine compute_projection_matrix


!-----------------------------------------------------------------------------!

  subroutine compute_model_perturbation(g,gt,cdi,cmi,dtd,gamma,a,epsilon,eta,dm,conf)

    use my_types
    use m_linalg

    implicit none

    ! subroutine arguments
    type(crs_matrix):: g,gt,cdi,cmi,dtd
    type(inv_conf):: conf
    type(vector) :: dm,gamma
    real(kind=dbl):: epsilon,eta
    type(matrix):: a

    ! local variables
    type(matrix):: at
    ! some temporary storage gets only allocated when needed
    type(crs_matrix) :: ctmp1
    type(crs_matrix) :: ctmpa,ctmpb,ctmpc,ctmpd
    type(matrix) :: mtmp1,mtmpa,mtmp2,mtmpb,mtmp3
    type(vector) :: vtmp1,vtmp2
    type(vector) :: vtmpa,vtmpb,vtmpc

    ! dm=-A*[a^T(Gcd^-1G+epsilon*Cm^-1+eta*D^t*D) A]^-1* A^t*gamma

    !G^T*Cd^-1 *G
    ctmp1=crsmatmul(gt,cdi)
    ctmpa=crsmatmul(ctmp1,g)
    call deallocate_crs_matrix(ctmp1)

    ! epsilon* Cm^-1
    ctmpb=crsscamatmul(epsilon,cmi)


    ! eta*D^T*D
    ctmpc=crsscamatmul(eta,dtd)

    ! Gcd^-1G+epsilon*Cm^-1+eta*D^t*D
    ctmp1=crsmatadd(ctmpa,ctmpb)
    ctmpd=crsmatadd(ctmp1,ctmpc)

    call deallocate_crs_matrix(ctmpa)
    call deallocate_crs_matrix(ctmpb)
    call deallocate_crs_matrix(ctmpc)
    call deallocate_crs_matrix(ctmp1)

    ! A^T*...*A
    mtmp1=crsmatmmatmul(ctmpd,a)
    at=mtranspose(a)
    mtmp2=mchangesignmat(at)
    mtmpa=mmatmul(at,mtmp1)

	call deallocate_crs_matrix(ctmpd)
    call deallocate_matrix(mtmp1)
    call deallocate_matrix(mtmp2)

    ! [...]^-1
    mtmpb=svd_inverse(mtmpa)
    call deallocate_matrix(mtmpa)

    ! -A* ...*A^T*gamma
    mtmp1=mchangesignmat(a)
    mtmp2=mmatmul(mtmp1,mtmpb)
    mtmp3=mmatmul(mtmp2,at)

    dm=mmatvecmul(mtmp3,gamma)

call deallocate_matrix(mtmpb)
call deallocate_matrix(mtmp1)
call deallocate_matrix(mtmp2)
call deallocate_matrix(mtmp3)



call deallocate_matrix(at)
  end subroutine compute_model_perturbation


!---------------------------------------------------------------------------------!

  subroutine use_pseudo_inverse(g,dobs,dpred,dm,conf)
    ! this subroutine uses svd for directly solving g*dm=dd
    use my_types
    use my_functions
    use m_linalg

    ! subroutine arguments
    type(crs_matrix):: g
    type(vector) :: dobs,dpred,dm
    type(inv_conf):: conf

    ! local varaibles
    integer ::i,j,k
    type(vector) :: vtmp1
    type(matrix) :: mtmp1,mtmp2,mtmp3
    type(matrix):: a,gi

    real(kind=dbl),allocatable:: u(:,:),s(:,:),v(:,:)

    type(matrix):: up,spi,vp

    integer :: nze

    a=crs2mat(g)

    call svd_decomposition(a%val,u,s,v)

nze=0
    do i=1,a%n2
       if (abs(s(i,i))>=conf%svdthre) then
          nze=nze+1
       end if
    end do

    if (nze<a%n2) then
       write(cdum,'(I3.3)') nze
       cdum= 'using the first '// cdum //' eigenvalues'
       call genmsg('m_inverse - pseudo_inverse',&
            'G is signular',cdum,2)
    end if



    up%n1=a%n1
    up%n2=nze
    allocate(up%val(up%n1,up%n2))
    up%val=0.0

    spi%n1=nze
    spi%n2=nze
    allocate(spi%val(spi%n1,spi%n2))
    spi%val=0.0

    vp%n1=a%n2
    vp%n2=nze
    allocate(vp%val(vp%n1,vp%n2))
    vp%val=0.0


    i=0
    do j=1,a%n2
       if (s(j,j)>=conf%svdthre) then
          i=i+1
          spi%val(i,i)=1.0/s(j,j)
          up%val(:,i)=u(:,j)
          vp%val(:,i)=v(:,j)
       end if
    end do

    ! pseudo inverse
    ! gi=v s^-1 u^t

    mtmp1=mmatmul(vp,spi)
    mtmp2=mtranspose(up)
    gi=mmatmul(mtmp1,mtmp2)

    call deallocate_matrix(mtmp1)
    call deallocate_matrix(mtmp2)

    ! dm=gi d

    vtmp1=vecsub(dobs,dpred)

    dm=mmatvecmul(gi,vtmp1)
call deallocate_matrix(gi)

  end subroutine use_pseudo_inverse

  !-------------------------------------------------------------------------------!

  subroutine update_velocity_model(vcur,vpert,dm)
    use my_types
    use my_functions
    implicit none

    ! subroutine arguments
    type(scalar_field_2d) :: vcur
    type(scalar_field_2d) :: vpert
    type(vector):: dm

    ! local variables
    integer :: i,j,id

    vpert%x0=vcur%x0;vpert%y0=vcur%y0
    vpert%nx=vcur%nx;vpert%ny=vcur%ny
    vpert%dx=vcur%dx;vpert%dy=vcur%dy
    vpert%cn=vcur%cn
    allocate(vpert%val(vpert%nx+2*vpert%cn,vpert%ny+2*vpert%cn))
    vpert%val=0.0_dbl

    do i=1,vcur%nx+vcur%cn*2
       do j=1,vcur%ny+vcur%cn*2
          call velnod2id(i,j,vcur,id)
          vpert%val(i,j)=dm%val(id)
          vcur%val(i,j)=vcur%val(i,j)+dm%val(id)
       end do
    end do

  end subroutine update_velocity_model

  !------------------------------------------------------------------------------!

  subroutine compute_raypath_coverage(g,vcur,rapco)

    use my_types
    use my_functions
    implicit none

    ! subroutine arguments
    type(crs_matrix)::g
    type(scalar_field_2d):: vcur,rapco

    ! local variables
    integer :: i,ii,jj

    rapco%x0=vcur%x0;rapco%y0=vcur%y0
    rapco%nx=vcur%nx;rapco%ny=vcur%ny
    rapco%dx=vcur%dx;rapco%dy=vcur%dy
    rapco%cn=vcur%cn

    allocate(rapco%val(rapco%nx+rapco%cn*2,rapco%ny+rapco%cn*2))
    rapco%val=0

    do i=1,size(g%col,1)
       call id2velnod(g%col(i),vcur,ii,jj)
       rapco%val(ii,jj)=rapco%val(ii,jj)+1
    end do

  end subroutine compute_raypath_coverage


end module m_inverse
