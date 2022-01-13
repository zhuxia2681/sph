module sli1

  use number
  use numerical

  contains

subroutine ferdi(l,m,wan,ngrid,offset,stepvec,Rnl)
!
!JDY: 2021.2.22
!call Ferdis relation.
!
! Rnl(rvec) = int(dOmega wnR(rvec) Y^*lm(theta,phi)
!
! where Omega is solid angle.
!

  implicit none

!  complex(8),allocatable,intent(in) :: wan(:,:,:)
  real(8),allocatable,intent(in) :: wan(:,:,:)
  real(8),intent(in) :: stepvec(3,3),offset(3)
  integer(4),intent(in) :: l,m,ngrid(3)

  complex(8),allocatable,intent(out) :: Rnl(:,:,:)

  integer(4) :: i,j,k,nx,ny,nz
  real(8) :: x,y,z,r,theta,phi
  complex(8) :: valZlm,valylm
  real(8),allocatable :: gridvec(:,:,:,:)
  
  allocate(Rnl(ngrid(1),ngrid(2),ngrid(3)))

  call grid_vec(ngrid,stepvec,offset,gridvec)

  do i = 1, ngrid(1)
    do j = 1, ngrid(2)
      do k = 1, ngrid(3)
         x = gridvec(i,j,k,1)
         y = gridvec(i,j,k,2)
         z = gridvec(i,j,k,3)
        call cart2polar(x,y,z,r,theta,phi) 
!        call ylm(l,m,theta,phi,valylm)
        call Zlm(l,m,theta,phi,valZlm)
!        Rnl(i,j,k) = wan(i,j,k)*valylm
        Rnl(i,j,k) = wan(i,j,k)*valZlm
      end do
    end do
  end do

  deallocate(gridvec)

end subroutine ferdi

subroutine grid_vec(ngrid,stepvec,offset,outvec)

  implicit none

  real(8),intent(in) :: stepvec(3,3),offset(3)
  real(8) :: x,y,z
  real(8),allocatable,intent(out) :: outvec(:,:,:,:)
  integer(4),intent(in) :: ngrid(3)
  integer(4) :: i,j,k,nx,ny,nz,m1,m2,m3,ng(3)

  ng = ngrid
!for even number grid, we must set -1 because points are 2n+1.
  if (mod(ng(1),2)==0) then
    m1 = 1
  else
    m1 = 0
  end if
  if (mod(ng(2),2)==0) then
    m2 = 1
  else
    m2 = 0
  end if
  if (mod(ng(3),2)==0) then
    m3 = 1
  else
    m3 = 0
  end if


  nx = int(ng(1)/2)
  ny = int(ng(2)/2)
  nz = int(ng(3)/2)


  allocate(outvec(ng(1),ng(2),ng(3),3))

  do i = -nx , nx - m1
    do j = -ny , ny - m2
      do k = -nz , nz - m2
        x = dble(i)*stepvec(1,1) + dble(i)*stepvec(2,1) + &
            dble(i)*stepvec(3,1) + offset(1)
        y = dble(j)*stepvec(1,2) + dble(j)*stepvec(2,2) + & 
            dble(j)*stepvec(3,2) + offset(2)
        z = dble(k)*stepvec(1,3) + dble(k)*stepvec(2,3) + &
            dble(k)*stepvec(3,3) + offset(3)
        outvec(i+nx+1,j+ny+1,k+nz+1,1) = x
        outvec(i+nx+1,j+ny+1,k+nz+1,2) = y
        outvec(i+nx+1,j+ny+1,k+nz+1,3) = z
      end do
    end do
  end do

end subroutine grid_vec

end module sli1
