subroutine rand(rnd)
 
  implicit none

  real(8),intent(out) :: rnd
  integer(4) :: i,seedsize
  integer(4),allocatable :: seed(:)
  
  call random_seed(size=seedsize)
  allocate(seed(seedsize))

  do i = 1,seedsize
    call system_clock(count=seed(i))
  end do
  
  call random_seed(put=seed(:))

  call random_number(rnd)

end subroutine rand
