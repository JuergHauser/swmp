! * m_inout *
!
! This module contains all the subroutines which are dealing with input and output operations
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

module m_inout

  implicit none

contains

  !---------------------------------------------------------------------------------!
  ! subroutines for reading .in files                                               !
  !---------------------------------------------------------------------------------!

  subroutine read_creobs_conf(conf)
    use my_types
    implicit none

    type(creobs_conf):: conf

    open(unit=input,file='creobs.in',status='old')
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%ifn_traveltimes
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%do_ran
    read(input,*) conf%stdev
    read(input,*) conf%rseed
    read(input,*) conf%unc
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%ofn_traveltimes
    close(input)
  end subroutine read_creobs_conf


!-------------------------------------------------------------------------------------!

  subroutine read_gcray_conf(ifn_gcray,sous,recs,conf)

    use my_types
    use my_constants
    implicit none
  character(len=*) :: ifn_gcray
    type(gcray_conf):: conf
    type(receivers):: recs
    type(sources):: sous

    integer :: i

    open(unit=input,file=ifn_gcray,status='old')

    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,'(a14)') conf%ifn_sources
    read(input,'(a14)') conf%ifn_receivers
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%vel
    read(input,*) rearth
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%ofn_traveltimes
    close(input)

    ! read sources
    open (unit=input,file=conf%ifn_sources,status='old')
    read(input,*) sous%nn
    allocate(sous%pos(sous%nn,2))
    do i=1,sous%nn
       read(input,*) sous%pos(i,1),sous%pos(i,2)
    end do
    close(input)


    ! read receivers
    open (unit=input,file=conf%ifn_receivers,status='old')
    read(input,*) recs%nn

    recs%mar=1
    allocate(recs%pos(recs%nn,2))
    allocate(recs%nar(recs%nn))
    allocate(recs%ariv(recs%nn,recs%mar))


    ! load receiver positions
    do i=1,recs%nn
       read(input,*) recs%pos(i,1),recs%pos(i,2)
    end do

    close(input)

  end subroutine read_gcray_conf

  !------------------------------------------------------------------------------------------------!

  subroutine read_gen2dv_conf(ifn_gend2dv,conf,spikes)
    use my_types
    implicit none
    character(len=*) :: ifn_gend2dv
    type(gen2dv_conf):: conf
    type(list_ind_val)::spikes
    integer :: i,ii
    character(len=strlen) :: line

    ! read in the configuration
    open(unit=input,file=ifn_gend2dv,status='old')
    read(input,*) cdum; read(input,*)  cdum; read(input,*) cdum
    read(input,*) conf%x0, conf%y0
    read(input,*) conf%nx, conf%ny
    read(input,*) conf%dx, conf%dy
    read(input,*) conf%cn
    read(input,*) cdum; read(input,*)  cdum; read(input,*) cdum
    read(input,*) conf%vbg
    read(input,*) cdum; read(input,*)  cdum; read(input,*) cdum
    read(input,*) conf%do_rand
    read(input,*) conf%stadev
    read(input,*) conf%rseed
    read(input,*) cdum;read(input,*) cdum; read(input,*) cdum
    read(input,*) spikes%n
    if (spikes%n>0) then
       allocate(spikes%ind(spikes%n,2))
       allocate(spikes%val(spikes%n))
       do i=1,spikes%n
          read(input,*) spikes%ind(i,1), spikes%ind(i,2), spikes%val(i)
       end do
    end if
    read(input,*) cdum;read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%do_cebo
    read(input,*) conf%cebo_pert
    read(input,*) conf%cebo_size
    read(input,*) conf%cebo_spac
    read(input,*) cdum;read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%do_modcov
    read(input,*) conf%covbg
    read(input,*) cdum;read(input,*) cdum; read(input,*) cdum
    read(input,'(a)') line
    conf%ofn_velmod = trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_modcov = trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_animod = trim(line(1:index(line,'!')-1))

    close(input)

  end subroutine read_gen2dv_conf

  !-----------

  subroutine read_inv_conf(ifn_inv,conf,epsilon,eta)
    ! this subrourine read sthe parameters for the inversion
    use my_types
    use my_constants
    implicit none

    ! subroutine arguments
    character(len=*) :: ifn_inv
    type(inv_conf) :: conf
    real(kind=dbl):: epsilon, eta
    ! local variables

    open(unit=input,file=ifn_inv,status='old')
    read(input,*) cdum; read(input,*) cdum;read(input,*) cdum
    read(input,*) conf%ifn_vinit
    read(input,*) conf%ifn_vcur
    read(input,*) conf%ifn_vunc
    read(input,*) conf%ifn_predtimes
    read(input,*) conf%ifn_obstimes
    read(input,*) conf%ifn_frechet_hdr
    read(input,*) conf%ifn_frechet_mat
    read(input,*) conf%ifn_frechet_rai
    read(input,*) cdum; read(input,*) cdum;read(input,*) cdum
    read(input,*) conf%invmod
    read(input,*) conf%subdim
    read(input,*) conf%svdthre
    read(input,*) epsilon
    read(input,*) conf%do_sds
    read(input,*) eta
    read(input,*) rearth
    read(input,*) cdum; read(input,*) cdum;read(input,*) cdum
    read(input,*) conf%ofn_vpert
    read(input,*) conf%ofn_vupd
    read(input,*) conf%ofn_tres
    read(input,*) conf%ofn_rapco
    close(input)

  end subroutine read_inv_conf

  !-------

  subroutine read_misfit_conf(ifn_misfit,conf)
    use my_types
    implicit none
    character(len=*) :: ifn_misfit
    type(misfit_conf):: conf
    open(unit=input,file=ifn_misfit,status='old')
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%ifn_vinv
    read(input,*) conf%ifn_vtrue
    read(input,*) conf%ifn_vref
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%nrx,conf%nry
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%ofn_dabs
    read(input,*) conf%ofn_drel
    read(input,*) conf%ofn_adrel
    read(input,*) conf%ofn_sumfile
    close (input)
  end subroutine read_misfit_conf

!---------------------------------------------------------------!

subroutine read_modro_conf(ifn_modro,conf)
use my_types
implicit none
 character(len=*) :: ifn_modro
type(modro_conf):: conf

open(unit=input,file=ifn_modro,status='old')

 read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
 read(input,*) conf%ifn_vcur
 read(input,*) conf%ifn_vref
 read(input,*) conf%nrx,conf%nry
 read(input,*) conf%rearth
end subroutine read_modro_conf

  !-------------------------------

subroutine read_rabgv_conf(ifn_rabgv,conf)

    use my_types
    implicit none

    character(len=*) :: ifn_rabgv
    type(rabgv_conf):: conf

open(unit=input,file=ifn_rabgv,status='old')
read(input,*) cdum;read(input,*) cdum;read(input,*) cdum
read(input,*) conf%ifn_mod
read(input,*) conf%ifn_rac
read(input,*) cdum;read(input,*) cdum;read(input,*) cdum
read(input,*) conf%bgv
read(input,*) conf%mnr
read(input,*) cdum;read(input,*) cdum;read(input,*) cdum
read(input,*) conf%ofn_mod
close(input)

end subroutine read_rabgv_conf

  !-------------------------------

  subroutine read_rat_conf(ifn_rat,sous,recs,recmode,wafc,vmod,conf)

    ! this subroutine reads in the parameters given in rat.in. It sets an intial
    ! value to the variables and allocate arrays given in the user derived types wafs,
    ! rex,vmod and params.

    use my_types
    use my_constants
    implicit none

    !subroutine arguments
    character(len=*)     :: ifn_rat
    type(sources)         :: sous
    type(wavefront_container):: wafc
    type(receivers)       :: recs
    type(receiver_mode):: recmode
    type(rat_conf)       :: conf
    type(scalar_field_2d) :: vmod

    ! local variables
    character(len=strlen)      :: line
    integer            :: i,j,k
    integer,pointer :: recsel(:)

    ! read rat.in
    open(unit=input,file=ifn_rat,status='old')
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
     read(input,'(a)') line
     conf%ifn_velmod = trim(line(1:index(line,'!')-1))
     read(input,'(a)') line
     conf%ifn_sources= trim(line(1:index(line,'!')-1))
     read(input,'(a)') line
     conf%ifn_receivers = trim(line(1:index(line,'!')-1))
    ! read parameters
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%wt%dt
    read(input,*) conf%wt%maxit
    read(input,*) conf%wt%n_senod
    read(input,*) conf%wt%mode
    read(input,*) conf%wt%solver
    read(input,*) conf%wt%interp
    read(input,*) recs%mar
    read(input,*) conf%velint
    read(input,*) rearth
    read(input,*) conf%recmode
    read(input,'(a)') line
    conf%ifn_recmode=trim(line(1:index(line,'!')-1))
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%out_wafint
    read(input,*) conf%do_rays
    read(input,*) conf%out_rays
    read(input,*) conf%do_frechet
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,'(a)') line
    conf%ofn_arrivals=trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_wavefronts=trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_raypaths=trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_frechet_hdr=trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_frechet_mat=trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_frechet_rai=trim(line(1:index(line,'!')-1))
    read(input,'(a)') line
    conf%ofn_sumfile=trim(line(1:index(line,'!')-1))
    close(input)

    ! read sources
    open (unit=input,file=conf%ifn_sources,status='old')
    read(input,*) sous%nn
    allocate(sous%pos(sous%nn,2))
    do i=1,sous%nn
       read(input,*) sous%pos(i,1),sous%pos(i,2)
    end do
    close(input)


    ! read receivers
    open (unit=input,file=conf%ifn_receivers,status='old')
    read(input,*) recs%nn

    allocate(recs%pos(recs%nn,2))
    allocate(recs%nar(recs%nn))
    recs%nar=0
    allocate(recs%stat(recs%nn))
    allocate(recs%ariv(recs%nn,recs%mar))

    ! load receiver positions
    do i=1,recs%nn
       read(input,*) recs%pos(i,1),recs%pos(i,2)
    end do
    close(input)
    recs%nar=0

	!print *,recs%nn


    wafc%n=conf%wt%maxit
    allocate(wafc%waf(conf%wt%maxit+1))
    wafc%waf(:)%n_nod=0

    ! load the velocity model

    call read_scalar_field_2d(vmod,conf%ifn_velmod)

    conf%xmin=vmod%x0-vmod%cn*vmod%dx
    conf%ymin=vmod%y0-vmod%cn*vmod%dy
    conf%xmax=(vmod%nx+vmod%cn)*vmod%dx+vmod%x0
    conf%ymax=(vmod%ny+vmod%cn)*vmod%dy+vmod%y0


    ! setup structure holding number of arrivals per receiver
    allocate(conf%tonar(recs%mar))
    conf%tonar=0

    if (conf%recmode==1) then

       allocate(recsel(recs%nn))

       ! read receiver selection

       open(unit=input,file=conf%ifn_recmode,status='old')

       recmode%nn=sous%nn
       allocate(recmode%act(recmode%nn))

       do i=1,recmode%nn
		  !print *,i
          read(input,*) recmode%act(i)%n,(recsel(k),k=1,recmode%act(i)%n)
          allocate( recmode%act(i)%val(recmode%act(i)%n))
          recmode%act(i)%val=recsel(1:recmode%act(i)%n)
       end do

       close(input)
       deallocate(recsel)
    end if


  end subroutine read_rat_conf

  !----------------------------------------------------------------------------------

  subroutine read_resa2d_conf(ifn_resa2d,conf)

    use my_types
    implicit none
   character(len=*) :: ifn_resa2d
    type(resa2d_conf):: conf


    ! read in the configuration
    open(unit=input,file=ifn_resa2d,status='old')
    read(input,*) cdum; read(input,*)  cdum; read(input,*) cdum
    read(input,*) conf%igrid
    read(input,*) cdum; read(input,*)  cdum; read(input,*) cdum
    read(input,*) conf%nrx, conf%nry
    read(input,*) cdum; read(input,*)  cdum; read(input,*) cdum
    read(input,'(a14)') conf%ovef
    read(input,'(a14)') conf%ogaf
    read(input,'(a14)') conf%gmtvef
    read(input,'(a14)') conf%gmtgaf
    read(input,*) conf%sw
    close(input)

  end subroutine read_resa2d_conf

  !--------------------------------------------------------------

  subroutine read_tresid_conf(ifn_tresid,conf)

    use my_types
    implicit none
    character(len=*) :: ifn_tresid
    type(tresid_conf):: conf

    open(unit=input,file=ifn_tresid,status='old')
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%ifn_predtimes
    read(input,*) conf%ifn_obstimes
    read(input,*) cdum; read(input,*) cdum; read(input,*) cdum
    read(input,*) conf%ofn_sumfile
    read(input,*) conf%dump_res
    read(input,*) conf%ofn_residuals

    close(input)

  end subroutine read_tresid_conf


!------------------------------------------------------------------
  subroutine read_vn2vj_conf(ifn_vn2vj,conf)
    use my_types
    implicit none

    !subroutine arguments
       character(len=*) :: ifn_vn2vj
    type(vn2vj_conf):: conf

    open(unit=input,file=ifn_vn2vj,status='old')
    read(input,*) cdum; read(input,*) cdum;read(input,*) cdum
    read(input,*) conf%ifn_vmod
    read(input,*) conf%cn
    read(input,*) cdum; read(input,*) cdum;read(input,*) cdum
    read(input,*) conf%ofn_vmod
    close(input)

  end subroutine read_vn2vj_conf


  !----------------------------------------------------------------------------------!
  ! subroutines for handling the frechet matrix
  !----------------------------------------------------------------------------!

  subroutine create_frechet_files(freha,matfile,raifile)
    ! This subroutine creates the binary and header file used for the
    ! storage of the Frechet Matrix.

    use my_types
    implicit none

    ! subroutine arguments
    ! type(rat_conf):: conf
    type(frechet_matrix_header):: freha
    character(len=*) matfile,raifile

    open(unit=output,file=matfile,status='replace',action='write',&
         access='direct',form='unformatted',recl=freha%recole)
    close(output)

    open(unit=output,file=raifile,status='replace')
    close(output)

  end subroutine create_frechet_files

  !----------------------------------------------------------------------------!

  subroutine append_frechet_row(frero,freha,matfile,raifile)
    ! This subroutine adds a row to the binary file and header file for the
    ! Frechet Matrix.

    use my_types
    implicit none

    ! subroutine arguments
    type(frechet_matrix_row):: frero
    type(frechet_matrix_header):: freha
    !type(rat_conf):: conf

    character(len=*) matfile,raifile

    ! local variables
    integer :: i
    integer:: nre

    nre=0

    open(unit=bin,file=matfile,status='old',action='write',&
         access='direct',form='unformatted',recl=freha%recole)
    do i=1,freha%n2
       if (frero%nze(i)==1) then
          nre=nre+1
          freha%nre=freha%nre+1
          write(bin,rec=freha%nre) freha%n1,i,frero%val(i)

       end if
    end do
    close (bin)


    open(unit=output,file=raifile,status='old',access='append')
    write(output,*)  freha%nre-nre+1, freha%nre
    close (output)


  end subroutine append_frechet_row

  !---------------------------------------------------------------------------!

  subroutine read_frechet_header(freha,hdrfile)
    ! write the header information for the frechet matrix file.
    use my_types
    implicit none

    ! subroutine arguments
    type(frechet_matrix_header):: freha

    character(len=*) :: hdrfile

    ! local varaibles
    integer :: i,n

    open(unit=input,file=hdrfile,status='old')
    read(input,*) freha%n1,freha%n2
    read(input,*) freha%nre
    read(input,*) freha%recole

    close(input)

  end subroutine read_frechet_header

  !----------------------------------------------------------------------------------!

  subroutine read_frechet_matrix(g,usarr,matfile,hdrfile,raifile)

    use my_types
    implicit none

    ! subroutine arguments
    type(crs_matrix) :: g
    type(list_index):: usarr
    ! local variables
    type(frechet_matrix_header):: freha

    integer :: ren1,ren2,nre
    integer:: i,j,k,icol,irow
    type(list_index):: rain

    type(matrix):: mtmp

    character(len=*) :: raifile,matfile,hdrfile

    ! read frechet matrix header infromation
    call read_frechet_header(freha,hdrfile)


    ! laod the ray index

    rain%n=freha%n1
    allocate (rain%ind(rain%n,2))
    open(unit=input,file=raifile,status='old')
    do i=1,rain%n
       read(input,*) rain%ind(i,1),rain%ind(i,2)
    end do
    close(input)

    ! count the number of elements needed for the frechet matrix.
    ! We do not now before hand which rows have to be loaded
    ! This is due to multiarrival tomography

    open(unit=input,file=raifile,status='old')

    j=1
    nre=0
    g%n1=0

    do i=1,rain%n
       if (usarr%ind(j,2)==i) then
          nre=nre+(rain%ind(i,2)-rain%ind(i,1))+1
          g%n1=g%n1+1
          j=j+1
       end if
    end do
    close(input)


    allocate(g%row(g%n1+1))
    allocate(g%val(nre))
    allocate(g%col(nre))


    ! open binary fileread(input,*) ren1,ren2
    open(unit=bin,file=matfile,status='old',action='read',&
         access='direct',form='unformatted',recl=freha%recole)

    icol=1
    irow=0
    j=1

    do i=1,rain%n
       if (usarr%ind(j,2)==i) then
          j=j+1

          irow=irow+1
          g%row(irow)=icol

          do k=rain%ind(i,1),rain%ind(i,2)
             read(bin,rec=k) idum,g%col(icol),g%val(icol)
             icol=icol+1
          end do

       end if
    end do

    close(bin)
    close(input)


    g%row(irow+1)=icol
    g%n2=freha%n2

!!$!!$
!!$    write(*,*) irow
!!$    pause
!!$
!!$    write(*,*) 'val'
!!$    write(*,*) g%val
!!$    write(*,*) 'col'
!!$    write(*,*) g%col
!!$write(*,*) g%col(59:61)
!!$    write(*,*) 'row'
!!$    write(*,*) g%row
!!$
!!$    write(*,*) 'size'
!!$    write(*,*) freha%n1,freha%n2

    call deallocate_list_index(rain)

  end subroutine read_frechet_matrix




  !---------------------------------------------------------------------------!

  subroutine write_frechet_header(freha,hdrfile)
    ! write the header information for the frechet matrix file.
    use my_types
    implicit none

    ! subroutine arguments
    type(frechet_matrix_header):: freha
    character(len=*) :: hdrfile

    open(unit=output,file=hdrfile,status='replace')
    write(output,*) freha%n1,freha%n2
    write(output,*) freha%nre
    write(output,*) freha%recole

  end subroutine write_frechet_header





  !
  ! subroutines for reading data
  !--------------------------------------------------------------------------------!
  subroutine read_observed_arrivals(obs,filename)
    ! this subroutine loads the traveltimes

    use my_types
    implicit none

    ! subroutine arguments
    type(observed_arrivals):: obs
    character(len=*) ::  filename

    ! local varaibles
    integer :: i

    open(unit=input,file=filename,status='old',action='read')
    read(input,*) obs%n
    ! allocate arrays
    allocate(obs%sou(obs%n))
    allocate(obs%rec(obs%n))
    allocate(obs%arn(obs%n))
    allocate(obs%art(obs%n))
    allocate(obs%azi(obs%n))

    allocate(obs%unc(obs%n))

    do i=1,obs%n
   !!  read(input,'(3i8,4f16.8)') obs%sou(i),obs%rec(i),obs%arn(i),obs%art(i),obs%azi(i),obs%unc(i)
       read(input,*) obs%sou(i),obs%rec(i),obs%arn(i),obs%art(i),obs%azi(i),obs%unc(i)

    end do

    close (input)

  end subroutine read_observed_arrivals

  !----------------------------------------------------------------------------------!

  subroutine read_predicted_arrivals(pred,filename)
    ! this subroutine loads the traveltimes

    use my_types
    implicit none

    ! subroutine arguments
    type(predicted_arrivals):: pred
    character(len=*) ::  filename

    ! local varaibles
    integer :: i

    open(unit=input,file=filename,status='old',action='read')
    read(input,*) pred%n
    ! allocate arrays
    allocate(pred%sou(pred%n))
    allocate(pred%rec(pred%n))
    allocate(pred%arn(pred%n))
    allocate(pred%art(pred%n))
    allocate(pred%azi(pred%n))
    do i=1,pred%n
       read(input,'(3i8,2f16.8)') pred%sou(i),pred%rec(i),pred%arn(i),pred%art(i),pred%azi(i)
    end do

    close (input)

  end subroutine read_predicted_arrivals




  !-----------------------------------------------------------------------------------!

  subroutine read_scalar_field_2d(vmod,filename)
    ! this subroutine load an intial velocity model
    use my_types
    implicit none
    ! subroutine arguments
    character(len=*) ::  filename
    type(scalar_field_2d) :: vmod

    ! local variables
    integer::i,j

    open(unit=input,file=filename,status='old',action='read')
    read(input,*) vmod%x0,vmod%y0
    read(input,*) vmod%nx,vmod%ny
    read(input,*) vmod%dx,vmod%dy
    read(input,*) vmod%cn
    allocate(vmod%val(vmod%nx+2*vmod%cn,vmod%ny+2*vmod%cn))
    do i=1,vmod%nx+vmod%cn*2
       do j=1,vmod%ny+vmod%cn*2
          read(input,*) vmod%val(i,j)
       end do
    end do
    close(input)

  end subroutine read_scalar_field_2d



  !--------------------------------------------------------------
  subroutine read_scalar_field_2d_vn(vmod,filename)
    ! this subroutine load an intial velocity model with flipped y
    use my_types
    implicit none
    ! subroutine arguments
    character(len=*) ::  filename
    type(scalar_field_2d) :: vmod

    ! local variables
    integer::i,j

    open(unit=input,file=filename,status='old',action='read')
    read(input,*) vmod%ny,vmod%nx
    read(input,*) vmod%y0,vmod%x0
    read(input,*) vmod%dy,vmod%dx
    ! read(input,*) vmod%cn
   ! vmod%cn=1
  !  vmod%nx=vmod%nx-2
  !  vmod%ny=vmod%ny-2
    allocate(vmod%val(vmod%nx+2*vmod%cn,vmod%ny+2*vmod%cn))
    ! do j=vmod%ny+vmod%cn*2,1,-1
    do i=1,vmod%nx+vmod%cn*2
       do j=vmod%ny+vmod%cn*2,1,-1
          read(input,*) vmod%val(i,j)
       end do
    end do
    close(input)

  end subroutine read_scalar_field_2d_vn



  !---------------------------------------------------------------------
  ! subroutines for writing data
  !---------------------------------------------------------------------

  subroutine write_arrival_data(recs,filename,sid,dor,spread)
    ! this subroutine writes the arrival times to a file

    use my_types
    implicit none

    ! subroutine arguments
    type(receivers),intent(in) :: recs
    character(len=strlen):: filename
    integer :: sid
integer :: dor,spread
    ! local variables
    integer :: h,i,j

    if (sid==1) then
       open(unit=output,file=filename,form='formatted',status='replace')
    else
       open(unit=output,file=filename,form='formatted',status='old',position='append')
    end if

if (spread==0) then
  do i=1,recs%nn
       do j=1,recs%nar(i)
          if( recs%ariv(i,j)%paic==1.or.dor==0) then
             write(output,'(i8,i8,i8,f16.8,f16.8,9f16.8)') sid,i,j,recs%ariv(i,j)%time,recs%ariv(i,j)%azi
          end if
       end do
    end do

else
    do i=1,recs%nn
       do j=1,recs%nar(i)
          if( recs%ariv(i,j)%paic==1.or.dor==0) then
             write(output,'(i8,i8,i8,f16.8,f16.8,9f16.8)') sid,i,j,recs%ariv(i,j)%time,recs%ariv(i,j)%azi,recs%ariv(i,j)%spf
          end if
       end do
    end do
end if


    close(output)
  end subroutine write_arrival_data

  !-------------------------------------------------------------------------!

  subroutine write_arrival_times(recs,filename,sid)
    ! this subroutine writes the arrival times to a file

    use my_types
    implicit none

    ! subroutine arguments
    type(receivers),intent(in) :: recs
    character(len=strlen):: filename
    integer :: sid
    ! local variables
    integer :: h,i,j

    if (sid==1) then
       open(unit=output,file=filename,form='formatted',status='replace')
    else
       open(unit=output,file=filename,form='formatted',status='old',position='append')
    end if
    do i=1,recs%nn
       do j=1,recs%nar(i)
          if( recs%ariv(i,j)%paic==1) then
             write(output,'(i8,i8,i8,f16.8)') sid,i,j,recs%ariv(i,j)%time
          end if
       end do
    end do

    close(output)
  end subroutine write_arrival_times

  !-------------------------------------------------------------------------!

subroutine write_observed_arrivals(obs,filename)

use my_types
implicit none

type(observed_arrivals):: obs
character(len=*) :: filename
integer :: i

open(unit=output,file=filename,status='replace')
write(output,*) obs%n
do i=1,obs%n
       write(output,'(3i8,4f16.8)') obs%sou(i),obs%rec(i),obs%arn(i),obs%art(i),obs%azi(i),obs%unc(i)
end do
close(output)

end subroutine write_observed_arrivals

!---------------------------------------------------------

subroutine write_predicted_arrivals(pred,filename)

use my_types
implicit none

type(predicted_arrivals):: pred
character(len=*) :: filename
integer :: i

open(unit=output,file=filename,status='replace')
write(output,*) pred%n
do i=1,pred%n
       write(output,'(3i8,9f16.8)') pred%sou(i),pred%rec(i),pred%arn(i),pred%art(i),pred%azi(i)
end do
close(output)

end subroutine write_predicted_arrivals






  !----------------------------------------------------------------------------!

subroutine write_misfit_summary(rms_abs,cor_abs,cor_ano,conf)

use my_types
implicit none

type(misfit_conf):: conf
real(kind=dbl)   :: rms_abs,cor_abs(3),cor_ano(3)

open(unit=output,file=conf%ofn_sumfile,status='replace')

write(output,*) 'model comaprison'
write(output,'(f12.6,A)') rms_abs, ' rms'
write(output,'(f12.6,A)') cor_abs(1), ' correlation coefficient'
write(output,'(f12.6,A)') cor_abs(2), ' significance level for rejecting H0'
write(output,'(f12.6,A)') cor_abs(3), ' fisher z'

write(Output,*) 'model perturbation comparison'
write(output,'(f12.6,A)') cor_ano(1), ' correlation coefficient'
write(output,'(f12.6,A)') cor_ano(2), ' significance level for rejecting H0'
write(output,'(f12.6,A)') cor_ano(3), ' fisher z'

close(output)

end subroutine write_misfit_summary

!-------------------------------------------------------------------------!

  subroutine write_tresid_summary(resid,conf)
    ! write information about the inversion to a file
    use my_types
    implicit none

    ! subroutine arguments
    type(tresid_conf):: conf
    type(traveltime_residuals) :: resid

    !local variables
       integer:: i

    open(unit=output,file=conf%ofn_sumfile,status='replace')

    write(output,'(f12.6,A)') resid%trms, '  rms travel time residuals'
  write(output,'(f12.6,A)') resid%tmean, '  mean travel time residuals'
  write(output,'(f12.6,A)') resid%tvar, '  variance travel time residuals'
    write(output,'(i12,A)') resid%tnar,'  number of arrivals'
    write(output,'(i12,A)') resid%n, ' set of arrivals'
    write(output,*) 'arn    nar    rms    mean    var'
    do i=1,resid%n
       write(output,'(i6,i6,f12.6,f12.6,f12.6)') i,int(resid%nar%val(i)),resid%rms%val(i),resid%mean%val(i), resid%var%val(i)
    end do
    close(output)

  end subroutine write_tresid_summary

  !-----------------------------------------------------------------------------!

  subroutine write_rat_summary(vmod,params,recs)
    ! this subroutine write the grid settign to a file which is than
    ! loaded by gmt for plotting the results
    use my_types
    implicit none

    ! subroutine arguments
    type(rat_conf):: params
    type(scalar_field_2d) :: vmod
    type(receivers) :: recs
    !local variables
    !    integer:: output
    integer ::i

    open(unit=output,file=params%ofn_sumfile,form='formatted')
    write(output,*) (vmod%nx+vmod%cn*2)*(vmod%ny+vmod%cn*2), '  number of velocity nodes'
    write(output,*) sum(params%tonar(:)) , '  number of raypaths'

    do i=1,recs%mar
       write(output,'(A,I3,A,i5)') ' arrival type ', i, '  number of arrivals' ,params%tonar(i)
    end do
    close(output)

  end subroutine write_rat_summary

  !---------------------------------------------------------------------------!

  subroutine write_raypaths(recs,conf)
    ! storing the raypaths to the disk. raypaths.dat contains all the raypaths while
    ! while raypaths01.dat, raypaths02.dat and raypaths03.dat contain
    ! the second and third arrivals

    use my_types
    implicit none

    ! subroutine arguments
    type(receivers) :: recs
    type(rat_conf) :: conf
    ! local variables
    integer :: i,j,k
    character(2) :: filenumber
    character(35) :: filename


    if (conf%out_rays==1.or.conf%out_rays==3) then

       if (conf%sid==1) then
          open(unit=output,file=conf%ofn_raypaths,form='formatted',status='replace')
       else
          open(unit=output,file=conf%ofn_raypaths,form='formatted',status='old',position='append')
       end if

       do i=1,recs%nn
          if (recs%nar(i)>0) then
             do j=1,recs%nar(i)
                if (recs%ariv(i,j)%paic==1) then
                   do k=1,recs%ariv(i,j)%path%n
                      write(output,*) recs%ariv(i,j)%path%pos(k,1),recs%ariv(i,j)%path%pos(k,2)
                   end do
                   write(output,'(a1)') '>'
                end if
             end do
          end if
       end do
       write(output,'(a1)') '>'
       close(output)

    end if

    if (conf%out_rays==2.or.conf%out_rays==3) then
       do j=1,recs%mar
          write(filenumber,'(I2.2)') j
          filename = conf%ofn_raypaths(1:index(conf%ofn_raypaths,'.')-1) &
          // '.' // filenumber // '.' &
          // conf%ofn_raypaths(index(conf%ofn_raypaths,'.')+1:index(conf%ofn_raypaths,' '))

          if (conf%sid==1) then
             open(unit=output,file=filename,form='formatted',status='replace')
          else
             open(unit=output,file=filename,form='formatted',status='old',position='append')
          end if
          do i=1,recs%nn
             if (recs%nar(i)>=j) then
                if (recs%ariv(i,j)%paic==1) then
                   do k=1,recs%ariv(i,j)%path%n
                      write(output,*) recs%ariv(i,j)%path%pos(k,1),recs%ariv(i,j)%path%pos(k,2)
                   end do
                   write(output,'(a1)') '>'
                end if
             end if
          end do
          !      write(output,'(a1)') '>'
          close(output)

       end do
    end if

  end subroutine write_raypaths
  !--------------------------------------------------------------------------------!

!!$  subroutine write_rays(wafs,conf)
!!$    ! writitng the wavefronts segmetns to a file The character > is
!!$    ! used for indicating that there is no connection between the segments
!!$    ! this allows easy plotting in gmt with psxy
!!$
!!$    use my_types
!!$    implicit none
!!$
!!$    ! subroutine arguments
!!$    type(wavefronts)     :: wafs
!!$    type(rat_conf)     :: conf
!!$    ! local variables
!!$    integer :: i,j,k
!!$    real(kind=dbl):: x,y
!!$
!!$
!!$    open(unit=output,file='ray.dat',form='formatted',status='replace')
!!$
!!$    do j=1,conf%wt%n_nod
!!$       do i=1,conf%wt%maxit
!!$          if (wafs%stat(j,i)==1) then
!!$             x=wafs%pos(j,i,1)
!!$             y=wafs%pos(j,i,2)
!!$
!!$             write(output,*) x,y
!!$          end if
!!$       end do
!!$       write(output,'(a1)') '>'
!!$    end do
!!$
!!$
!!$    close(output)
!!$  end subroutine write_rays

  !----------------------------------------------------------------------------------!

  subroutine write_residuals(resid,conf)

    use my_types
    implicit none
! subroutine arguments
    type(traveltime_residuals):: resid
    type(tresid_conf):: conf

! local variables
    character(2) :: filenumber
    character(35) :: filename
integer :: i,j


    if (conf%dump_res==1.or.conf%dump_res==3) then
       filename=conf%ofn_residuals
       open (unit=output,file=conf%ofn_residuals,status='replace')

       do i=1,resid%n
          do j=1,resid%ares(i)%n
             write(output,*) resid%ares(i)%val(j)
          end do
       end do
       close(output)
    end if

    if (conf%dump_res==1.or.conf%dump_res==3) then
       do i=1,resid%n
          write(filenumber,'(I2.2)') i
          filename= 'arn'//filenumber//conf%ofn_residuals
          open(unit=output,file=filename,status='replace')
          do j=1,resid%ares(i)%n
             write(output,*) resid%ares(i)%val(j)
          end do
          close(output)
       end do

    end if

  end subroutine write_residuals

!------------------------------------------------------------------------------!

  subroutine write_sia_summary(arms,brms,drms,p,siadec,conf)
    use my_types
    implicit none

    type(sia_conf):: conf
    real(kind=dbl)::arms,brms,drms,p
    integer :: siadec

    open(unit=output,file=conf%ofn_sumfile,status='replace')

    write(output,'(f12.6,A)') conf%tk,' temperature'
    write(output,'(f12.6,A)') brms,' rms model b'
    write(output,'(f12.6,A)') arms,' rms model a'
    write(output,'(f12.6,A)') drms,' rms difference'
    if (siadec>1) then
       write(output,'(f12.6,A)') p, ' p'
       write(output,'(f12.6,A)') exp(-drms/conf%tk), ' rho'

    end if
    if (siadec==1) then
       write(output,*) 'uncondtional acceptance'
       write(output,*) 'brms < arms'
    else if (siadec==2) then
       write(output,*) 'conditional acceptance'
       write(output,*) 'brms > arms & p < rho'
    else
       write(output,*) 'rejection'
       write(output,*) 'brms > arms & p >= rho'
    end if

  end subroutine write_sia_summary

!-----------------------------------------------------------------------------!

  subroutine write_scalar_field_2d(vmod,filename)
    ! this subroutine load an intial velocity model
    use my_types
    implicit none
    ! subroutine arguments
    character(len=*) ::  filename
    type(scalar_field_2d) :: vmod

    ! local variables
    integer::i,j

    open(unit=output,file=filename,status='replace')
    write(output,*) vmod%x0,vmod%y0
    write(output,*) vmod%nx,vmod%ny
    write(output,*) vmod%dx,vmod%dy
    write(output,*) vmod%cn
    do i=1,vmod%nx+vmod%cn*2
       do j=1,vmod%ny+vmod%cn*2
          write(output,*) vmod%val(i,j)
       end do
    end do
    close(output)
  end subroutine write_scalar_field_2d

!--------------------------------------------------------------------------

subroutine write_gmt_scalar_field_2d(vmod,filename,headername)
use my_types
implicit none
    character(len=*) ::  filename
    character(len=*) :: headername
    type(scalar_field_2d) :: vmod
integer :: i,j

    open(unit=output,file=filename,status='unknown')
     do i=1+vmod%cn,vmod%nx+vmod%cn
        do j=1+vmod%cn,vmod%ny+vmod%cn
           write(output,*) vmod%val(i,j)
        end do
     end do
close(output)

     open(unit=output,file=headername,form='formatted')
     write(output,*) vmod%x0
     write(output,*) vmod%x0+(vmod%nx-1)*vmod%dx
     write(output,*) vmod%y0
     write(output,*) vmod%y0+(vmod%ny-1)*vmod%dy
     write(output,*) vmod%dx
     write(output,*) vmod%dy
     close(output)



end subroutine write_gmt_scalar_field_2d

  !-----------------------------------------------------------------------------------!

subroutine write_sia_tmp(conf)
use my_types
implicit none
type(sia_conf):: conf

open (unit=output,file='sia.tmp',status='replace')
write(output,*) conf%k
write(output,*) conf%tk
write(output,*) conf%lk
write(output,*) conf%ak
write(output,*) conf%u
write(output,*) conf%spert
write(output,*) conf%sacep
close(output)

end subroutine write_sia_tmp


!-------------------------

  subroutine write_vector(vec,filename)
    use my_types
    implicit none
    ! subroutine arguments
    character(len=*) ::  filename
    type(vector) :: vec

    ! local variables
    integer :: i


    open(unit=output,file=filename,status='replace')
    write(output,*) vec%n
    do i=1,vec%n
       write(output,*) vec%val(i)
    end do

  end subroutine write_vector

  !-----------------------------------------------------------------------------------!


  subroutine write_wavefronts(wafc,conf)
    ! writitng the wavefronts segmetns to a file The character > is
    ! used for indicating that there is no connection between the segments
    ! this allows easy plotting in gmt with psxy

    use my_types
    use my_functions
    implicit none

    ! subroutine arguments
   ! type(wavefronts)     :: wafs
    type(wavefront_container):: wafc
    type(rat_conf)     :: conf
    type(node):: cur
    ! local variables
    integer :: i,j,k,h,m
    ! real(kind=dbl):: x,y


    if (conf%sid==1) then
       open(unit=output,file=conf%ofn_wavefronts,form='formatted',status='replace')
    else
       open(unit=output,file=conf%ofn_wavefronts,form='formatted',status='old',position='append')
    end if

    do i=1,conf%wt%maxit
       if (mod(i,conf%out_wafint)==0.and.wafc%waf(i)%n_nod>0) then

          k=wafc%waf(i)%headpid
          j=1




          do while (j<=wafc%waf(i)%n_nod)
       !   write(*,*) k,h
             call locate_in_list_index(wafc%waf(i)%index,k,h)
                   !write(*,*)      wafc%waf(i)%index%ind(h,2)
             write(output,*) wafc%waf(i)%pos(wafc%waf(i)%index%ind(h,2),1:2)
             if (wafc%waf(i)%con(wafc%waf(i)%index%ind(h,2),2)==0) then
                write(output,'(a1)') '>'
             end if
             k=wafc%waf(i)%neig(wafc%waf(i)%index%ind(h,2),2)
             j=j+1
          end do

          k=wafc%waf(i)%headpid
          call locate_in_list_index(wafc%waf(i)%index,k,h)
          write(output,*) wafc%waf(i)%pos(wafc%waf(i)%index%ind(h,2),1:2)

          write(output,'(a1)') '>'
       end if
    end do
    write(output,'(a1)') '>'
    close(output)
  end subroutine write_wavefronts

end module m_inout
