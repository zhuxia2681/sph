module numerical

  use number

  implicit none

  contains

subroutine Zlm(l,m,theta,phi,val)
!Real spherical Harmonics
!JDY: 2021.2.22

  implicit none

  real(8),intent(in) :: theta,phi
  integer(4),intent(in) :: l,m

  complex(8),intent(out) :: val
  complex(8)  :: tmp1,tmp2
  

  if (m == 0) then
    call ylm(l,0,theta,phi,val)
  elseif (m > 0) then
    call ylm(l,-m,theta,phi,tmp1)
    call ylm(l,m,theta,phi,tmp2)
    val = (1/sqrt(2.0)) * (tmp1 + ((-1)**m) * &
           tmp2)
  elseif (m < 0) then
    call ylm(l,m,theta,phi,tmp1)
    call ylm(l,-m,theta,phi,tmp2)
    val = (IMG/sqrt(2.0)) * (tmp1 - ((-1)**m) * &
           tmp2)
  end if

end subroutine Zlm


subroutine ylm(l,m,theta,phi,val)

  implicit none
  
  real(8),intent(in) :: theta,phi
  integer(4),intent(in) :: l,m

  integer(4) :: tmpnume, tmpdeno
  real(8)  :: sgn, coeff, nume, deno, test, tmp

  complex(8),intent(out) :: val

  if (m .gt. 0) then
    sgn = (-1)**m
  else if (m .le. 0) then
    sgn = 1
  end if

!call Legendre polynomial
  call plm(l,m,theta,tmp)

  call fact(l-abs(m),tmpnume)
  call fact(l+abs(m),tmpdeno)

  nume = tmpnume
  deno = tmpdeno

  coeff = sqrt((2 * l + 1)/(4 * pi)) * sqrt((nume)/(deno))
  val = sgn * coeff * exp(IMG * m * phi) * tmp

end subroutine ylm

subroutine plm(l,m,theta,val)

  implicit none

  integer(4),intent(in) ::  l, m
  real(8),intent(in) :: theta
  real(8),intent(out) :: val

  if (l .eq. 0 .and. m .eq. 0) then
    val = 1.0
  else if (l .eq. 1 .and. abs(m) .eq. 1) then
    val = sin(theta)
  else if (l .eq. 1 .and. m .eq. 0) then
    val = cos(theta)
  else if (l .eq. 2 .and. abs(m) .eq. 2) then
    val = 3 * sin(theta) ** 2
  else if (l .eq. 2 .and. abs(m) .eq. 1) then
    val = 3 * sin(theta) * cos(theta)
  else if (l .eq. 2 .and. m .eq. 0) then
    val = 0.5* (3 * sin(theta) ** 2 - 1)
  end if

end subroutine plm



subroutine fctrl(x,f)
!
! calculate factorial of x
!
  implicit none

  real(8),intent(inout) :: x
  real(8),intent(out) :: f 
!local
  integer(4) :: i,imax

  f = 1.0
  if (abs(x) .le. 0.1) x = 0.0
  if (x .lt. 0) then
    f = 0.0
  elseif (x .gt. 0) then
    imax = x + 0.1
    do i = 1, imax
      f = f * real(i)
    end do
  end if
  return

end subroutine fctrl

subroutine fact(n, x)

  implicit none
  integer(4), intent(in) :: n
  integer(4) :: i
  integer(4), intent(out) :: x

  x = 1
  if (n == 0) then
    x = 1
  else
    do i=1, n
      x = x * i
    end do
  end if

end subroutine fact

subroutine cart2polar(x,y,z,r,theta,phi)
 
  implicit none

  real(8),intent(in) :: x,y,z
  real(8),intent(out) :: r,theta,phi
  real(8) :: r2

  r2 = x**2 + y**2 + z**2
  r = sqrt(r2)
  if (r < EPS9) then
    theta = PI / 2
  else
    theta = acos(z / r)
  end if

  if (x > EPS9) then
    phi = atan(y/x)
  elseif (x < -1 * EPS9) then
    phi = atan(y/x) * PI
  else
    phi = sin(y) * (PI/2.0)
  end if

end subroutine cart2polar

subroutine polar2cart(r,theta,phi,x,y,z)

  implicit none

  real(8),intent(in) :: r,theta,phi
  real(8),intent(out) :: x,y,z

  x = r * sin(theta) * cos(phi)
  y = r * sin(theta) * sin(phi)
  z = r * cos(theta)

end subroutine polar2cart

subroutine conv2bohr(sfo,stepvec,sfoout,stepvecout)

  implicit none

  real(8),intent(in) :: sfo(3),stepvec(3,3)
  real(8),intent(out) :: sfoout(3),stepvecout(3,3)
!JDY local
  integer(4) :: i,j

  do i = 1, 3
    sfoout(i) = sfo(i) / BOHR
  end do

  do i = 1,3
    do j = 1,3
      stepvecout(i,j) = stepvec(i,j) / BOHR
    end do
  end do

end subroutine conv2bohr


subroutine double_to_int(dd,ii)

  implicit none

  real(8),intent(in) :: dd
  integer(4),intent(out) :: ii

  if (dd < 0.0) then
    ii = int(dd - 0.5)
  else
    ii = int(dd + 0.5)
  end if

end subroutine double_to_int



end module numerical

