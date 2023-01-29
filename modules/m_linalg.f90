! * m_linalg *
!
! contains function for linear algebra operations using matrices stored
! in as arrays or in the crs notation. 
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

module m_linalg

  implicit none

contains

  !---------------------------------------------------------------------------------!

  function mtranspose(a) result(b)
    ! wrap around the intrinsic tranpose function
    use my_types
    implicit none
    type(matrix) :: a,b
    b%n1=a%n2
    b%n2=a%n1

    allocate(b%val(b%n1,b%n2))
    b%val=transpose(a%val)

  end function mtranspose

  !----------------------------------------------------------------------------------!

  function mmatmul(a,b) result(c)
    ! wrap around the intrinsic matmul function
    use my_types
    implicit none
    type(matrix) :: a,b,c
    c%n1=a%n1
    c%n2=b%n2
    allocate (c%val(c%n1,c%n2))
    c%val=matmul(a%val,b%val)

  end function mmatmul

  !------------------------------------------------------------------------------------!

  function mmatvecmul(a,b) result(c)
    ! wrap around the intrinsic matmul function
    use my_types
    implicit none
    type(matrix) :: a
    type(vector) :: b,c
    c%n=a%n1
    allocate (c%val(c%n))
    c%val=matmul(a%val,b%val)

  end function mmatvecmul

  !--------------------------------------------------------------------------------!

  function mvecmatmul(a,b) result(c)

    use my_types
    implicit none
    type(matrix) :: b
    type(vector) :: a,c

    c%n=b%n2
    allocate (c%val(c%n))
    c%val=matmul(a%val,b%val)

  end function mvecmatmul

  !-----------------------------------------------------------------------------------!

  function mchangesignmat(a) result(b)
    ! changes the sign of the matrix
    use my_types
    implicit none

    ! subroutine arguments
    type(matrix) :: a,b

    b%n1=a%n1
    b%n2=a%n2

    allocate(b%val(b%n1,b%n2))
    b%val=-a%val

  end function mchangesignmat

  !-------------------------------------------------------------------------------------!

  function mat2crs(mat,thre) result(crs)

    ! this subroutine converts a sparse matrix into a matrix stored
    ! in the compact row storage format only elements > thre are 
    ! retained

    use my_types
    implicit none

    ! subroutine arguments
    type(matrix)     :: mat
    type(crs_matrix) :: crs
    real(kind=dbl):: thre

    ! local variables
    integer :: i,j,s,t
    real(kind=dbl) :: a(mat%n1,mat%n2),val(mat%n1*mat%n2)
    integer :: row(mat%n1),col(mat%n1*mat%n2)

    s=0 ! column counter
    t=0 ! row counter

    a=mat%val

    do i=1,mat%n1       ! loop over columns
       if (count(abs(a(i,:))>thre)==0) then
          write(*,*) 'ERR: mat2crs no nonzero elements in a column'
          write(*,*) ' storage in crs format not possible'
          stop
       end if
       do j=1,mat%n2
          if (abs(a(i,j))>thre) then
             s=s+1
             val(s)=a(i,j)
             col(s)=j
             if (i/=t) then
                t=t+1
                row(t)=s
             end if
          end if
       end do
    end do


    ! allocate the crs matrix

    allocate(crs%val(s))
    allocate(crs%col(s))
    allocate(crs%row(t+1))

    crs%val(:)=val(1:s)
    crs%col(:)=col(1:s)
    crs%row(1:t)=row(1:t)
    crs%row(t+1)=s+1
    crs%n1=mat%n1
    crs%n2=mat%n2


  end function mat2crs

  !---------------------------------------------------------------------------------!

  function crs2mat(crs) result(mat)
    ! converts a crs sotred matrix into a standard notation
    use my_types
    implicit none

    ! subroutine arguments
    type(matrix)     :: mat
    type(crs_matrix) :: crs

    ! local variables
    integer :: i,j

    allocate(mat%val(crs%n1,crs%n2))
    mat%val=0.0_dbl

    j=0
    do i=1,size(crs%val,1)
       if (j+1<=size(crs%row)) then      ! evt -1
          if (crs%row(j+1)==i) then
             j=j+1    
          end if
       end if
       mat%val(j,crs%col(i))=crs%val(i)
    end do

    mat%n1=crs%n1
    mat%n2=crs%n2

  end function crs2mat

  !-----------------------------------------------------------------------------!



  function crstranspose(a) result(b)
    ! this subroutien computes the transpose for a matrix in crs storage
    ! subroutine arguments
    use my_types
    implicit none

    ! subroutine arguments
    type(crs_matrix):: a,b

    integer :: colcount(a%n2)

    integer i,j,k,k1

    allocate(b%val(size(a%val)))
    allocate(b%col(size(a%col)))
    allocate(b%row(a%n2+1))


    ! local variables

    ! count number of elemts n each column of A
    colcount=0
    do k=1,size(a%val)
       colcount(a%col(k))=colcount(a%col(k))+1
    end do

    ! knowing the number of element in each colum of A allows
    ! the construction ofr b%row

    b%row(1)=1
    do j=1,a%n2
       b%row(j+1)=b%row(j)+colcount(j)
    end do

    ! we have setup b%row and next assign b%col and b%val, by running thorught 
    ! b%val in the natular order for a (which is scatterte order for B).
    ! for j in 1:n colcount(j) will accumalte the number of elemts visited in 
    ! column j
    b%col=0
    b%val=0
    colcount=0

    do i=1,a%n1
       do k=a%row(i),a%row(i+1)-1
          j=a%col(k)
          ! compute the index k1 for A's k-th element in the array B
          k1=b%row(j)+colcount(j)
          ! assing b%val and b%col for this element
          b%val(k1)=a%val(k)  
          b%col(k1)=i
          ! update colcount
          colcount(j)=colcount(j)+1
       end do
    end do

    b%n1=a%n2
    b%n2=a%n1

  end function crstranspose

  !----------------------------------------------------------------------------!

  function crsinversdiag(a) result(b)
    ! tis function computes the inverse of a square diagonal crs matrix. This is straight
    ! forward as the invers of a diagonal matrix is given by one over its diagonal
    ! elements.
    use my_types
    implicit none

    ! subroutine arguments
    type(crs_matrix) :: a,b
    ! local variables
    integer :: i

    allocate(b%val(size(a%val)))
    allocate(b%col(size(b%val)))
    allocate(b%row(size(a%row)))

    b%n1=a%n2
    b%n2=a%n1

    b%col=a%col
    b%row=a%row
    do i=1,size(a%val)
       b%val=1.0_dbl/a%val(i)
    end do

  end function crsinversdiag

  !------------------------------------------------------------------------------!

  function crsmatmul(a,b) result(d)
    ! this subroutine performs the matrix multpilication for two sparse matrices 
    !
    ! for i=1 to n       loop over rows of A
    ! for j=1 to n       loop over cols of A 
    ! for k=1 to n       loop over cols of B 
    ! c(i,k)=c(i,k)+a(i,j)*b(j,k)
    ! end for
    ! end for
    ! end for
    !
    ! An efficient method for performing the matrix-matrix multiplication, , 
    ! involved three nested loops, which essentially, find matching non-empty rows, 
    ! C(i,:) and A(i,:), find the matching row B(k,:) for each non-zero A(i,k), 
    ! and then find matching non-zero entries, C(i,j) and B(k,j).
    ! This can only be done efficiently if the matrix formats for A, B, and C 
    ! provide efficient access to non-empty rows

    use my_types
    implicit none

    ! subroutine arguments
    type(crs_matrix) :: a,b,d

    ! local variables
    type(crs_matrix) :: c
    integer        :: i,j,k,s,ss
    real(kind=dbl) :: temp(a%n1*a%n2)
    integer        :: ptr(a%n1*a%n2),itemp,ind(a%n1*a%n2),mlast,bigint,nlast

    ! first check if we can do the matrix multiplication if not
    ! we will run into a segmentation fault. The execution is therfore stopped.
    if (a%n2/=b%n1) then
       write(*,*) 'ERR: crsmatmul inner matrix dimensions must agree'
       stop
    end if

    bigint=a%n1*a%n2
    allocate(c%val(a%n1*b%n2))
    allocate(c%row(a%n1*b%n2))
    allocate(c%col(a%n1*b%n2))

    mlast=0
    nlast=0
    ! loop over the non zero rows of A
    do i=1,size(a%row,1)-1
       ptr=0
       itemp=0
       ind=0
       ! loop over the nonzero colum elements of the nonzero row in a
       do j=a%row(i),a%row(i+1)-1
          !loop over the colums of B
          do k=b%row(a%col(j)),b%row(a%col(j)+1)-1
             if( (ptr(b%col(k)) == 0))then
                itemp=itemp+1
                ptr(b%col(k))=itemp
                temp(ptr(b%col(k))) = a%val(j)*b%val(k)
                ind(ptr(b%col(k))) = b%col(k)
             else 
                temp(ptr(b%col(k))) =  temp(ptr(b%col(k)))+a%val(j)*b%val(k)
             end if

          end do
       end do
       !we need to write these results into a new spare matrix in crs storage notation
       !we have everything for the current row ;-)
       do s = 1,itemp
          ss=minloc(ind(1:itemp),1)  
          ! might be not so nice, but sorts columwise which needs to be the case for the successfull
          ! multiplication of the resulting crs matrix with another csr storage matrix.The resulting
          !  crs storage matrix is the formulated according to the definition of the crs storage scheme     
          c%val(mlast+s) = temp(ss)
          c%col(mlast+s) = ind(ss)
          ind(ss)=bigint            
       end do
       mlast = mlast + itemp
       nlast=nlast+1
       c%row(nlast)=mlast-itemp+1
    end do

    ! rewritte the results properly into sparse matrix, has to be done
    ! as we don't know the size of the result beforehand.
    allocate(d%val(mlast))
    allocate(d%col(mlast))
    allocate(d%row(nlast+1))

    d%val(:)=c%val(1:mlast)
    d%col(:)=c%col(1:mlast)
    d%row(1:nlast)=c%row(1:nlast)
    d%row(nlast+1)=mlast+1
    d%n1=a%n1
    d%n2=b%n2

call deallocate_crs_matrix(c)

  end function crsmatmul

  !----------------------------------------------------------------------------!

  function crsmatvecmul (a,x) result(y)
    ! computes the matrix vector product for a crs matrix witha a vector
    use my_types
    implicit none
    !subroutine arguments
    type(crs_matrix) :: a
    type(vector) :: x,y
    ! local variables
    integer :: i,j

    ! first check if we can do this operation form a dimension point of view
    if (a%n2/=x%n) then
       write(*,*) 'ERR: crsmatvecmul inner matrix dimensions must agree'
       stop
    end if

    y%n=a%n1
    allocate(y%val(y%n))
    y%val=0

    do i=1,a%n1
       do j=a%row(i),a%row(i+1)-1
          y%val(i)=y%val(i)+a%val(j)*x%val(a%col(j))
       end do
    end do

  end  function crsmatvecmul

  !------------------------------------------------------------------------------!


  function crsmatmmatmul(a,b) result(c)
    ! this function multiplies a crs format storage matrix wth a standard matrix.
    use my_types
    implicit none
    ! subroutine arguments
    type(crs_matrix) :: a
    type(matrix) :: b
    type(matrix) ::c

    !local variables
    integer        :: i,j,k
    integer :: jj

    ! for i=1 to n       loop over rows of A
    ! for j=1 to n       loop over cols of A and rwos of B
    ! for k=1 to n       loop over cols of B 
    ! c(i,k)=c(i,k)+a(i,j)*b(j,k)
    ! end for
    ! end for
    ! end for


    ! first check if we can do the amtrix multiplication if not
    ! we will run into a segmentation fault. The execution is therfore stopped.
    if (a%n2/=b%n1) then
       write(*,*) 'ERR: crsmatmul inner matrix dimensions must agree'
       stop
    end if

    !write(*,*) a%n1,a%n2
    !write(*,*) b%n1,b%n2
    !pause

    c%n1=a%n1
    c%n2=b%n2
    allocate(c%val(c%n1,c%n2))

    do i=1,size(a%row,1)-1
       ! loop over the nonzero colum elements of the nonzero row in a
       do j=a%row(i),a%row(i+1)-1
          ! loop over the corresponding elements in b     
          do k=1,b%n2
             jj=a%col(j)
             !    write(*,*) i,jj,k
             c%val(i,k)=c%val(i,k)+a%val(j)*b%val(jj,k)
          end do
       end do
    end do


  end function crsmatmmatmul

  !------------------------------------------------------------------------------!

  function crsmatadd(a,b) result(d)
    ! adds two crs matrixes
    use my_types
    implicit none

    ! subroutien arguments
    type(crs_matrix) :: a,b,d
    !local varibles
    type(crs_matrix) :: c
    integer :: i,ka,kb,kc,kamax,kbmax
    integer :: ja,jb

    allocate(c%val(size(a%val)+size(b%val)))
    allocate(c%col(size(a%col)+size(b%col)))
    allocate(c%row(size(a%row)))   


    ! first verify if we can add them

    if (a%n1/=b%n1.or.a%n2/=b%n2) then
       write(*,*) 'ERR: crsmatadd matrix dimension dont agree' 
       stop
    end if

    kc=1
    c%row(1)=kc

    ! loop over all the row in a and b
    do i=1,a%n1 
       ka=a%row(i)
       kb=b%row(i)
       kamax=a%row(i+1)-1
       kbmax=b%row(i+1)-1

       ! as fortran does not know repeat until we construct a repeat until structure
       ! using do exit (good old pascal would provide do until ;-) )

       do

          if ( ka <= kamax ) then
             ja = a%col(ka)
          else
             ja = a%n2 + 1
          end if

          if ( kb <= kbmax ) then
             jb = b%col(kb)
          else
             jb = b%n2 + 1
          end if


          if ( kamax < ka .and. kbmax < kb ) then
             exit ! finished this row
          end if

          !  Three cases

          if ( ja == jb ) then
             c%val(kc) = a%val(ka) +b%val(kb)
             c%col(kc) = ja
             ka = ka + 1
             kb = kb + 1
             kc = kc + 1
          else if ( ja < jb ) then
             c%col(kc) = ja
             c%val(kc) = a%val(ka)
             ka = ka + 1
             kc = kc + 1
          else if ( jb < ja ) then
             c%col(kc) = jb
             c%val(kc) = b%val(kb)
             kb = kb + 1
             kc = kc + 1
          end if

       end do

       c%row(i+1) = kc

    end do


    ! rewritte the results properly into sparse matrix, has to be done
    ! as we don't know the size of the result beforehand.
    allocate(d%val(kc-1))
    allocate(d%col(kc-1))
    allocate(d%row(a%n1+1))

    d%val(:)=c%val(1:kc-1)
    d%col(:)=c%col(1:kc-1)
    d%row(1:a%n1)=c%row(1:a%n1)
    d%row(a%n1+1)=a%n1+1
    d%n1=a%n1
    d%n2=b%n2
    
    call deallocate_crs_matrix(c)

  end function crsmatadd

  !------------------------------------------------------------------------------!
  function scavecmul(q,a) result(b)
    ! multipies a vector with a scalar or vice versa
    use my_types
    implicit none
    !subroutine arguments
    type(vector) :: a,b
    real(kind=dbl) :: q

    b%n=a%n
    allocate(b%val(b%n))

    b%val=a%val*q

  end function scavecmul

  !--------------------------------------------------------------------------------!

  function crsscamatmul(q,a) result(b)
    ! multiplies a scalar with a crs matrix
    use my_types
    implicit none
    ! subroutine arguments
    type(crs_matrix) :: a,b
    real(kind=dbl) :: q

    b%n1=a%n1
    b%n2=a%n2
    allocate(b%val(size(a%val)))
    allocate(b%col(size(a%col)))
    allocate(b%row(size(a%row)))

    b%val=a%val*q
    b%col=a%col
    b%row=a%row

  end function crsscamatmul

  !------------------------------------------------------------------------------!

  function vecsub(a,b) result(c)
    ! subtracts two vectors
    use my_types
    implicit none
    ! subroutine arguments
    type(vector) :: a,b
    type(vector) :: c

    if (a%n/=b%n) then
       write(*,*) 'ERR: vecsub vector dimensions must agree'
       stop
    end if

    c%n=a%n
    allocate(c%val(c%n))
    c%val=a%val-b%val

  end function vecsub

  !------------------------------------------------------------------------------

  function vecadd(a,b) result(c)
    ! adds two vectors
    use my_types
    implicit none
    ! subroutine arguments
    type(vector) :: a,b
    type(vector) :: c

    if (a%n/=b%n) then
       write(*,*) 'ERR: vecadd vector dimensions must agree'
       stop
    end if

    c%n=a%n
    allocate(c%val(c%n))
    c%val=a%val+b%val

  end function vecadd

  !----------------------------------------------------------------------------!

  function matadd(a,b) result(c)
    ! adds two matrices
    use my_types
    implicit none
    ! subroutine arguments
    type(matrix) :: a,b
    type(matrix) :: c

    if (a%n1/=b%n1.or.a%n2/=b%n2) then
       write(*,*) 'ERR: matadd matrix dimensions must agree'
       stop
    end if
    c%n1=a%n1
    c%n2=a%n2
    allocate(c%val(c%n1,c%n2))
    c%val=a%val+b%val

  end function matadd


  !----------------------------------------------------------------------------!
  function svd_orthogonalization (a,thre) result(b)
    ! othogonalization for a matrix using singular value decomposition. 
    ! the vectors are also normailzed.
    use my_types
    use my_functions
    implicit none

    ! subroutine arguments
    type(matrix) :: a,b
    real(kind=dbl):: thre
    ! local variabkes
    real(kind=dbl), allocatable:: u(:,:),s(:,:),v(:,:)
    integer :: i,j,k


!write(*,*) a%val

    ! decompose a into u,s,v
    call svd_decomposition(a%val,u,s,v)


    if (a%n2==1) then
       b%n1=a%n1
       b%n2=a%n2
       allocate(b%val(b%n1,b%n2))
       b%val=a%val

       call normalize_vector(a%val(:,1))

    else

       ! count vectors with wiegthing above threshold
       b%n2=0
       do i=1,a%n2
       !    write(*,*) i,s(i,i)
          if (abs(s(i,i))>=thre) then
             b%n2=b%n2+1
          end if
       end do

!!$ open(unit=10,file='b',status='replace')
!!$  do i=1,a%n1
!!$     write(10,'(50(32F32.16))') a%val(i,:)
!!$  end do
!!$  close(10)
       !write(*,*) 'sdim',b%n2

       if (b%n2<a%n2) then
          write(*,*) 'Due to redundancy, the subspace dimension'
          write(*,*) 'is being reduced from ',a%n2,' to ',b%n2
       end if

       b%n1=a%n1
       allocate(b%val(b%n1,b%n2))
       i=0
       do j=1,a%n2
          if (s(j,j)>=thre) then
             !  write(*,*) 'j',j
             i=i+1
             b%val(:,i)=u(:,j)
             
call normalize_vector(b%val(:,i))
          end if
       end do
    end if

!!if (b%n2<2) then
!!pause
!!end if
!!$u(:,1)=b%val(:,1)
!!$b%val(:,1)=b%val(:,2)
!!$b%val(:,2)=u(:,1)

    ! clean up
    deallocate(u,s,v)

  end function svd_orthogonalization


  !----------------------------------------------------------------------------!

  function  svd_inverse(a) result(ai)
    ! compute the inverse of a amtrix a using svd a=u*s*v^t and a^-1=v*s^-1*u^T
    ! crs_a and crs_ai are given ie. returned in crs storage. In the subroutien itself they are
    ! are covnerted into standard matrixes ie. cionverted into crs amtrices
    use my_types
    implicit none
    ! subroutien arguments
    type(matrix) :: a,ai
    ! local variables
    real(kind=dbl), allocatable:: u(:,:),s(:,:),v(:,:),si(:,:)
    integer :: i

    ! check for singularity 

    ! decompose a into u,s,v
    call svd_decomposition(a%val,u,s,v)


    ! we compute the inverse of the diagonal matrix s
    allocate(si(size(s,1),size(s,2)))
    si=0.0_dbl

    do i=1,size(s,1)
       si(i,i)=1.0_dbl/s(i,i)
    end do


    ! we build the matrix ai 

    ai%n1=a%n1
    ai%n2=a%n2

    allocate(ai%val(ai%n1,ai%n2))

ai%val=0.0_dbl

    ai%val=matmul(matmul(v,si),transpose(u))


  end function  svd_inverse

!-----------------------------------------------------------------------------!

  function lu_determinant(a) result(det)


    use my_types
    implicit none

    ! subroutine arguments
    type(matrix):: a
    real(kind=dbl):: det

    ! locla variables
    real(kind=dbl),allocatable:: lu(:,:)
    real(kind=dbl):: d
    integer :: i

    call lu_decomposition(a%val,lu,d)
    det=d

    do i=1,size(lu,1)
       det=det*lu(i,i)
    end do

  end function lu_determinant


  !------------------------------------------------------------------------------!

  subroutine svd_decomposition(a,u,s,v)
    ! This is just a wrap aorund svdcmp. It modifies the calling in such a way that a is not
    ! overwritten by u
    use my_types
    implicit none

    ! subroutine arguments
    real(kind=dbl) :: a(:,:)
    real(kind=dbl), allocatable :: u(:,:),s(:,:),v(:,:)

    ! local variables
    real(kind=dbl), allocatable :: w(:)
    integer :: i,m,n

    m=size(a,1);n=size(a,2)

    allocate(u(m,n))
u=0.0_dbl
    allocate(s(n,n))
s=0.0_dbl
    u=a

    call svdcmp(u,w,v)
    s=0.0_dbl
    do i=1,n
       s(i,i)=w(i)
    end do

  end subroutine svd_decomposition


!---------------------------------------------------------------------------!


  SUBROUTINE svdcmp(a,w,v)
    ! given a matrix a(1:m,1:n), this subroutine computes it singular value decomposition a=u*w*v^T. 
    ! The matrix u(1:m,1:n) replaces a on output. The diagonal matrix v is given as a vector. 
    ! The matrix v has the dimensions (n,n)

    use my_types
    use my_functions
    IMPLICIT NONE
    REAL(kind=dbl), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(kind=dbl), DIMENSION(:), INTENT(OUT),allocatable :: w
    REAL(kind=dbl), DIMENSION(:,:), INTENT(OUT),allocatable :: v
    INTEGER :: i,its,j,k,l,m,n,nm
    REAL(kind=dbl) :: anorm,c,f,g,h,s,scale,x,y,z
    REAL(kind=dbl), DIMENSION(size(a,1)) :: tempm
    REAL(kind=dbl), DIMENSION(size(a,2)) :: rv1,tempn
    m=size(a,1)
    n=size(a,2)
    allocate(v(n,n))
    allocate(w(n))
    g=0.0
    scale=0.0
    do i=1,n
       l=i+1
       rv1(i)=scale*g
       g=0.0
       scale=0.0
       if (i <= m) then
          scale=sum(abs(a(i:m,i)))
          if (scale /= 0.0) then
             a(i:m,i)=a(i:m,i)/scale
             s=dot_product(a(i:m,i),a(i:m,i))
             f=a(i,i)
             g=-sign(sqrt(s),f)
             h=f*g-s
             a(i,i)=f-g
             tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
             a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
             a(i:m,i)=scale*a(i:m,i)
          end if
       end if
       w(i)=scale*g
       g=0.0
       scale=0.0
       if ((i <= m) .and. (i /= n)) then
          scale=sum(abs(a(i,l:n)))
          if (scale /= 0.0) then
             a(i,l:n)=a(i,l:n)/scale
             s=dot_product(a(i,l:n),a(i,l:n))
             f=a(i,l)
             g=-sign(sqrt(s),f)
             h=f*g-s
             a(i,l)=f-g
             rv1(l:n)=a(i,l:n)/h
             tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
             a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
             a(i,l:n)=scale*a(i,l:n)
          end if
       end if
    end do
    anorm=maxval(abs(w)+abs(rv1))
    do i=n,1,-1
       if (i < n) then
          if (g /= 0.0) then
             v(l:n,i)=(a(i,l:n)/a(i,l))/g
             tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
             v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
          end if
          v(i,l:n)=0.0
          v(l:n,i)=0.0
       end if
       v(i,i)=1.0
       g=rv1(i)
       l=i
    end do
    do i=min(m,n),1,-1
       l=i+1
       g=w(i)
       a(i,l:n)=0.0
       if (g /= 0.0) then
          g=1.0_dbl/g
          tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
          a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
          a(i:m,i)=a(i:m,i)*g
       else
          a(i:m,i)=0.0
       end if
       a(i,i)=a(i,i)+1.0_dbl
    end do
    do k=n,1,-1
       do its=1,30
          do l=k,1,-1
             nm=l-1
             if ((abs(rv1(l))+anorm) == anorm) exit
             if ((abs(w(nm))+anorm) == anorm) then
                c=0.0
                s=1.0
                do i=l,k
                   f=s*rv1(i)
                   rv1(i)=c*rv1(i)
                   if ((abs(f)+anorm) == anorm) exit
                   g=w(i)
                   h=pythag(f,g)
                   w(i)=h
                   h=1.0_dbl/h
                   c= (g*h)
                   s=-(f*h)
                   tempm(1:m)=a(1:m,nm)
                   a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                   a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                end do
                exit
             end if
          end do
          z=w(k)
          if (l == k) then
             if (z < 0.0) then
                w(k)=-z
                v(1:n,k)=-v(1:n,k)
             end if
             exit
          end if


          if (its == 30) then
             call nrerror('m_linalg - svdcmp','no convergence')
          end if


          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dbl*h*y)
          g=pythag(f,1.0_dbl)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do j=l,nm
             i=j+1
             g=rv1(i)
             y=w(i)
             h=s*g
             g=c*g
             z=pythag(f,h)
             rv1(j)=z
             c=f/z
             s=h/z
             f= (x*c)+(g*s)
             g=-(x*s)+(g*c)
             h=y*s
             y=y*c
             tempn(1:n)=v(1:n,j)
             v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
             v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
             z=pythag(f,h)
             w(j)=z
             if (z /= 0.0) then
                z=1.0_dbl/z
                c=f*z
                s=h*z
             end if
             f= (c*g)+(s*y)
             x=-(s*g)+(c*y)
             tempm(1:m)=a(1:m,j)
             a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
             a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
          end do
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
       end do
    end do


  contains

    FUNCTION outerprod(a,b) result(op)
      use my_types
      REAL(kind=dbl), DIMENSION(:), INTENT(IN) :: a,b
      REAL(kind=dbl), DIMENSION(size(a),size(b)) :: op
      op = spread(a,dim=2,ncopies=size(b)) * &
           spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod


    real(kind=dbl) function pythag(a,b)
      use my_types
      implicit none

      !  real(kind=dbl) :: pythag
      real(kind=dbl) :: a,b
      real(kind=dbl) :: absa,absb

      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
         pythag=absa*sqrt(1.0_dbl+(absb/absa)**2)
      else
         if(absb.eq.0.0)then
            pythag=0.0
         else
            pythag=absb*sqrt(1.0_dbl+(absa/absb)**2)
         endif
      endif
      return
    end function pythag

  END SUBROUTINE svdcmp

!_-------------------------------------------------------------------------!

  subroutine lu_decomposition(a,lu,d)

    use my_types
    ! subroutine arguments
    real(kind=dbl):: a(:,:)
    real(kind=dbl),allocatable:: lu(:,:)
    real(kind=dbl):: d

    ! local varaibles
    integer :: n
    integer,allocatable:: indx(:)

    n=size(a,1)
    allocate(lu(n,n))
    lu=a
    call ludcmp(lu,indx,d)

  end subroutine lu_decomposition



	SUBROUTINE ludcmp(a,indx,d)
use my_types
use my_functions
	IMPLICIT NONE
	REAL(kind=dbl), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER, DIMENSION(:), allocatable :: indx
	REAL(kind=dbl), INTENT(OUT) :: d
	REAL(kind=dbl), DIMENSION(size(a,1)) :: vv
	REAL(kind=dbl), PARAMETER :: TINY=1.0e-20_sgl
	INTEGER :: j,n,imax
	n=size(a,1)

allocate(indx(n))

	d=1.0
	vv=maxval(abs(a),dim=2)

	if (any(vv == 0.0)) then
  call nrerror('m_linalg - svdcmp','singular matrix')
end if
	vv=1.0_sgl/vv
	do j=1,n
		imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
		if (j /= imax) then
			call swap(a(imax,:),a(j,:))
			d=-d
			vv(imax)=vv(j)
		end if
		indx(j)=imax
		if (a(j,j) == 0.0) a(j,j)=TINY
		a(j+1:n,j)=a(j+1:n,j)/a(j,j)
		a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
	end do

contains


  FUNCTION imaxloc(arr) result(res)
    use my_types
    REAL(kind=dbl), DIMENSION(:), INTENT(IN) :: arr
    INTEGER :: imax(1)
    integer :: res
    imax=maxloc(arr(:))
    res=imax(1)
  END FUNCTION imaxloc

     FUNCTION outerprod(a,b) result(op)
      use my_types
      REAL(kind=dbl), DIMENSION(:), INTENT(IN) :: a,b
    REAL(kind=dbl), DIMENSION(size(a),size(b)) :: op
      op = spread(a,dim=2,ncopies=size(b)) * &
           spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod


  subroutine swap(a,b)
      use my_types
      implicit none
      real(kind=dbl) :: a(:),b(:),dum(size(a))
      dum=a
      a=b
      b=dum
    end subroutine swap


	END SUBROUTINE ludcmp
!----------------------------------------------------------------------------------!




end module m_linalg
