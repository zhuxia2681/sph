program test

  implicit none

  real(8) :: rr
  integer(4) :: i
 
  do i = 1, 50
 
    call rand(rr)
    write(*,*) rr-0.5

  end do

end program test
