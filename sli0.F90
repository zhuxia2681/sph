subroutine sli0(l1,l2,ngrid,pnl1,pnl2,r,val)

  use number
  
  implicit none

  real(8),allocatable,intent(in) :: r(:,:,:)
  real(8),allocatable,intent(in) :: pnl1(:,:,:),pnl2(:,:,:)
  real(8),intent(out) :: val

  real(8),allocatable :: rkp1(:,:,:)
  real(8),allocatable :: f1(:,:,:),f2(:,:,:)
  integer(4),intent(in) :: ngrid(3),l1,l2
  integer(4) :: n1,ix,iy,iz,i,j,k,m,ng(3)
  integer(4) :: nokmax,ni
  real(8) :: a0,a1,a2,b0,b1,b2
  real(8) :: ho12,eras

  ng = ngrid

  do n1 = 1, 2

    do ix = 1, ng(1)
      do iy = 1, ng(1)
        do iz = 1, ng(1)
          f1(ix,iy,iz) = pnl1(ix,iy,iz) ** 2
          if (n1 <= 1) then 
            f2(ix,iy,iz) = f1(ix,iy,iz)
            k = abs(la - lb) - 2
          else
            f2(ix,iy,iz) = pnl2(ix,iy,iz) ** 2
            k = -2
          endif
        end do
      end do
    end do

    nokmax = min(la,lb) + 1
    allocate(rkp1(ngrid(1),ngrid(2),ngrid(3)))

    do n = 1, nokmax
      rkp1(1,1,1) = 0.001
      a1 = 0.0
      b1 = 0.0
      k = abs(k) + 2       
      do ix = 1, ng(1)
        do iy = 1, ng(1)
          do iz = 1, ng(1)
            if (n <= 1) then
              rkp1(ix,iy,iz) = r(ix,iy,iz)
              if (k <= 0) exit
                do m = 1, k
                  rkp1(ix,iy,iz) = rkp1(ix,iy,iz) * r(ix,iy,iz)
                end do
              end if
            else
              rkp1(ix,iy,iz) = r(ix,iy,iz) * r(ix,iy,iz)**2
            end if
          end do
        end do
      end do
!     |
! calc r1 integral
!     |
      xi(1) = 0.0
      xj(1) = 0.0
      xa(1) = 0.0
      xb(1) = 0.0
      a2 = 0.0    
      b2 = 0.0
!JDY assume all 2nd neighbor of origin 
!    is same distance.
      ho12 = r(1,1,2) / 12.0     

      ni = 40

      do ix = 3, ng(1), 2
        do iy = 3, ng(1), 2
          do iz = 3, ng(1), 2
            a0 = a2
            b0 = b2
            eras = 8.0 * f1(i)
          end do
        end do
      end do

end subroutine sli0
