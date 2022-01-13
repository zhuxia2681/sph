module number

  implicit none

  integer(4),parameter :: NCELL = 2
  integer(4),parameter :: NRAND = 2500000
  integer(4),parameter :: RAND_MAX = 2147483647
  real(8),parameter :: HARTREE = 27.21138505 
  real(8),parameter :: EPS9 = 1.0E-9
  real(8),parameter :: INF9 = 1.0E+9
  real(8),parameter :: BOHR = 0.52917721092
  real(8), parameter :: PI = 3.141592653589793238462643383279
  real(8), parameter :: FPI = 4 * PI
  complex(8),parameter :: IMG = (0.0,1.0)

end module number
