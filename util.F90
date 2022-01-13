module util

  implicit none
 
  contains

subroutine err(name,ierr)

  implicit none

  integer(4) :: ierr
  character(30) :: name
 
  if (ierr /= 0) then
    write(*,*) " Error at ",name," err code ",ierr
  endif
  
end subroutine err


subroutine clock_start(t)

  implicit none

  integer(4),intent(inout) :: t

  call system_clock(t)

end subroutine clock_start


subroutine clock_end(t1,t)

  implicit none

  integer(4),intent(in) :: t1
  real(8),intent(out) :: t
  integer(4) :: t_rate,t_max
  integer(4) :: diff,t2

  call system_clock(t2,t_rate,t_max)

  if (t2 < t1) then
    diff = (t_max - t1) + t2 + 1
  else
    diff = t2 - t1
  end if

  t = diff / dble(t_rate)
  
end subroutine clock_end

end module util
