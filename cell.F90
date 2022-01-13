module cell

  use number

  implicit none

contains

subroutine volume(vec,vol)

  implicit none

  real(8),intent(in) :: vec(3,3)
  real(8),intent(out) :: vol

  vol = abs(vec(1,1)*(vec(2,2)*vec(3,3)-vec(2,2)*vec(1,3))- &
            vec(1,2)*(vec(1,1)*vec(2,3)-vec(2,1)*vec(1,3))+ &
            vec(1,3)*(vec(1,1)*vec(2,2)-vec(2,1)*vec(1,2)))

end subroutine volume

subroutine cell_set(ngrid,uc,ucu,sfo,sfv,uco,ucv,&
           cellvec)
!
! Perform unit conversions and convert to cartesian 
! from lattice coordinates if necessary
!
!JDY:
!    sfo is scalar field origin
!    sfv is scalar field vector
!
! ucu is unit cell unit.
! ucu = 0 : bohr
!     = 1 : angstrom
!     = 2 : latvec

  implicit none

  logical,intent(in) :: uc
  real(8),intent(in) :: uco(3),ucv(3,3),sfo(3),sfv(3,3)
  real(8),intent(out) :: cellvec(3,3)
  integer(4),intent(in) :: ucu,ngrid(3)

!JDY local
  real(8) :: cellorigin(3)
  real(8) :: c1,c2,c3
  integer(4) :: i,j,k

  if (uc) then
! the case bohr
    if (ucu == 0) then
      cellorigin(1) = uco(1) * BOHR
      cellorigin(2) = uco(2) * BOHR
      cellorigin(3) = uco(3) * BOHR
      do i = 1, 3
        cellvec(i,1) = ucv(i,1) * BOHR
        cellvec(i,2) = ucv(i,2) * BOHR
        cellvec(i,3) = ucv(i,3) * BOHR
      end do
!the case latvec
    elseif (ucu == 2) then
      c1 = uco(1)
      c2 = uco(2)
      c3 = uco(3)
      cellorigin(1) = c1 * sfv(1,1) + c2 * sfv(2,1) + c3 * sfv(3,1)    
      cellorigin(2) = c1 * sfv(1,2) + c2 * sfv(2,2) + c3 * sfv(3,2)    
      cellorigin(3) = c1 * sfv(1,3) + c2 * sfv(2,3) + c3 * sfv(3,3)    
      do i = 1, 3
        c1 = ucv(i,1)
        c2 = ucv(i,2)
        c3 = ucv(i,3)
        cellvec(i,1) = c1 * sfv(1,1) + c2 * sfv(2,1) + c3 * sfv(3,1)
        cellvec(i,2) = c1 * sfv(1,2) + c2 * sfv(2,2) + c3 * sfv(3,2)
        cellvec(i,3) = c1 * sfv(1,3) + c2 * sfv(2,3) + c3 * sfv(3,3)
      end do
    end if
  else
    do i = 1, 3
      cellorigin(i) = sfo(i)
    end do
    do i = 1, 3
      do j = 1, 3
        cellvec(i,j) = sfv(i,j)
      end do
    end do
  end if

end subroutine cell_set

subroutine isovalue_scale(power,threshold,isovalue,sff,ngrid)

  implicit none

  real(8),allocatable,intent(in) :: sff(:,:,:)
  integer(4),intent(in) :: power,ngrid(3)
  real(8),intent(in) :: threshold
  real(8),intent(out) :: isovalue

  integer(4),parameter :: ns = 64
  real(8) :: samp,smin,smax,sless

! local
  integer(4) :: i,j,k,p,ss
  real(8) :: sa(ns),sb(ns),sc(ns),sd,se,sf
  real(8) :: w
!debug
  integer(4) :: ii
 
  smin = INF9
  sless = INF9
  smax = -INF9

  do i = 1, ngrid(1)
    do j = 1, ngrid(2)
      do k = 1, ngrid(3)
        if (sff(i,j,k) < smin) then
          smin = sff(i,j,k)
        end if
        if (sff(i,j,k) > smax) then
          smax = sff(i,j,k)
        end if
        if (abs(sff(i,j,k)) < sless) then
          sless = sff(i,j,k)
        end if
      end do  
    end do  
  end do  

!debug
  write(*,*) 
  write(*,*) " scalar value info. "
  write(*,*) " sless ",sless
  write(*,*) " smin ",smin
  write(*,*) " smax ",smax
  write(*,*) 
  
  if (abs(smax) > abs(smin)) then 
    samp = abs(smax)
  else
    samp = abs(smin)
  end if

  if (power > 0) then
    sa(1) = 0.0
    do ss = 1, ns-1
      sa(ss+1) = dble(ss)/dble(ns-1)
    end do
    do ss = 1, ns
      sb(ss) = sa(ss) * samp
    end do
    sb(1) = sb(1) - EPS9
    sb(ns) = sb(ns) + EPS9
    do ss = 1, ns
      sc(ss) = 0.0
    end do
    do i = 1, ngrid(1) 
      do j = 1, ngrid(2) 
        do k = 1, ngrid(3) 
          sd = abs(sff(i,j,k))
          se = 1.0
          do p = 1, power
            se = se * sd
          end do
          do ss = 1, ns
            if (sb(ss) < sd) then
              sc(ss) = sc(ss) + se
            end if
          end do
        end do
      end do
    end do

    sf = sc(1)

    if (sf > EPS9) then
      do ss = 1, ns
        sc(ss) = sc(ss) / sf
      end do

      call inversion(ns,sa,sc,threshold,w)

      isovalue = w
    else
      isovalue = 0.0
    end if
  end if

  isovalue = isovalue * samp

end subroutine isovalue_scale

subroutine inversion(n,xx,yy,zz,w)

  implicit none

  integer(4),intent(in) :: n
  real(8),intent(in) :: xx(n),yy(n),zz
  real(8),intent(out) :: w

!JDY
  integer(4) :: i,j
  real(8) :: a,b,c,s,t,v,u

  u = INF9
  j = 0

  do i = 1, n - 2
    v = (yy(i) - zz)**2 + (yy(i) - zz)**2 + (yy(i) - zz)**2
    if (v < u) then
      j = i
      u = v
    end if
  end do


  a = ((xx(j+1) - xx(j)) * (yy(j+2)-yy(j+1)) - &
       (xx(j+2) - xx(j+1)) * (yy(j+1) - yy(j)))/ &
       ((xx(j+2) - xx(j)) * (xx(j+2) - xx(j+1)) * (xx(j+1) - xx(j)))

  b = ((xx(j+1) - xx(j)) * (yy(j+2) - yy(j+1)) + &
       (xx(j+2) - xx(j+1)) * (yy(j+1) - yy(j))) / &
       (2.0 * (xx(j+2) - xx(j+1)) * (xx(j+1) - xx(j))) - &
       a * (xx(j) + 2.0 * xx(j+1) + xx(j+2)) / 2.0

  c = (yy(j) + yy(j+1) + yy(j+2 ) - b * (xx(j) + xx(j+1) + xx(j+2)) - &
      a * (xx(j)**2 + xx(j+1)**2 + xx(j+2)**2)) / 3.0


  u = (-b + sqrt(b**2 - 4.0 * a * (c - zz))) / (2.0 * a)
  v = (-b - sqrt(b**2 - 4.0 * a * (c - zz))) / (2.0 * a)

  s = (xx(j) - u)**2 + (xx(j+1) - u)**2 + (xx(j+2) - u)**2
  t = (xx(j) - v)**2 + (xx(j+1) - v)**2 + (xx(j+2) - v)**2


  if (s < t) then
    w = u
  else
    w = v
  end if

end subroutine inversion

subroutine trimming(ng,valin,ngout,valout)

  implicit none

  integer(4),intent(in) :: ng(3)
  integer(4),intent(out) :: ngout(3)
  real(8),allocatable,intent(in) :: valin(:,:,:)
  real(8),allocatable,intent(out) :: valout(:,:,:)
  integer(4) :: i,j,k

  allocate(valout(ng(1)-1,ng(2)-1,ng(3)-1))

  do i = 1, ng(1) -1
    do j = 1, ng(2) -1
      do k = 1, ng(3) - 1
        valout(i,j,k) = valin(i,j,k)
      end do  
    end do  
  end do  

  do i = 1, 3
    ngout(i) = ng(i) - 1
  end do

end subroutine trimming

end module cell
