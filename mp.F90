module mp

  implicit none

#ifdef PARA
include "mpif.h"
#endif

  public :: mp_bcast
 
  interface mp_bcast
    module procedure mp_bcast_int
    module procedure mp_bcast_real
    module procedure mp_bcast_double
    module procedure mp_bcast_cplx
    module procedure mp_bcast_char
    module procedure mp_bcast_logical
  end interface mp_bcast

  contains

subroutine mp_ini(ierr)

  implicit none

  integer(4),intent(out) :: ierr
  
  call MPI_init(ierr)
  if (ierr /= 0) then
    write(*,*) " error: mpi ini. failed. "
    return 
  endif

end subroutine mp_ini


subroutine mp_fin(ierr)

  implicit none

  integer(4),intent(out) :: ierr
  
  call MPI_finalize(ierr)
  if (ierr /= 0) then
    write(*,*) " error: mpi fin. failed. "
    return 
  endif

end subroutine mp_fin

subroutine mp_info(rank,size,ierr)

  implicit none

  integer(4),intent(out) :: rank,size,ierr

#ifdef PARA
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  if (ierr /= 0) then
    write(*,*) " error: mpi comm rank failed in mp_info. "
    return 
  endif
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  if (ierr /= 0) then
    write(*,*) " error: mpi comm size failed in mp_info. "
    return 
  endif
#else
  rank = 0
  size = 1
  ierr = 0
#endif

end subroutine mp_info

subroutine mp_bcast_int(val,n,rank,ierr)

  implicit none

  integer(4),intent(inout) :: val
  integer(4),intent(in) :: n,rank
  integer(4),intent(out) :: ierr

#ifdef PARA

  call MPI_Bcast(val,n,MPI_INTEGER,rank,MPI_COMM_WORLD,ierr)
 
#endif  

  return 

end subroutine mp_bcast_int

subroutine mp_bcast_real(val,n,rank,ierr)

  implicit none

  real(8),intent(inout) :: val
  integer(4),intent(in) :: n,rank
  integer(4),intent(out) :: ierr

#ifdef PARA

  call MPI_Bcast(val,n,MPI_REAL,rank,MPI_COMM_WORLD,ierr)
 
#endif  

  return 

end subroutine mp_bcast_real

subroutine mp_bcast_double(val,n,rank,ierr)

  implicit none

  real(16),intent(inout) :: val
  integer(4),intent(in) :: n,rank
  integer(4),intent(out) :: ierr

#ifdef PARA

  call MPI_Bcast(val,n,MPI_DOUBLE_PRECISION,rank,MPI_COMM_WORLD,ierr)
 
#endif  

  return 

end subroutine mp_bcast_double

subroutine mp_bcast_cplx(val,n,rank,ierr)

  implicit none

  complex(16),intent(inout) :: val
  integer(4),intent(in) :: n,rank
  integer(4),intent(out) :: ierr

#ifdef PARA

  call MPI_Bcast(val,n,MPI_COMPLEX,rank,MPI_COMM_WORLD,ierr)
 
#endif  

  return 

end subroutine mp_bcast_cplx

subroutine mp_bcast_logical(val,n,rank,ierr)

  implicit none

  logical,intent(inout) :: val
  integer(4),intent(in) :: n,rank
  integer(4),intent(out) :: ierr

#ifdef PARA

  call MPI_Bcast(val,n,MPI_LOGICAL,rank,MPI_COMM_WORLD,ierr)
 
#endif  

  return 

end subroutine mp_bcast_logical

subroutine mp_bcast_char(val,n,rank,ierr)

  implicit none

  character(256),intent(inout) :: val
  integer(4),intent(in) :: n,rank
  integer(4),intent(out) :: ierr

#ifdef PARA

  call MPI_Bcast(val,n,MPI_CHARACTER,rank,MPI_COMM_WORLD,ierr)
 
#endif  

  return 

end subroutine mp_bcast_char


subroutine mp_send(val,n,rank,tag,ierr)

  implicit none

  real(8),intent(inout) :: val
  integer(4),intent(in) :: rank,n,tag
  integer(4),intent(out) :: ierr
   

#ifdef PARA
  call MPI_Send(val,n,MPI_REAL,rank,tag,MPI_COMM_WORLD,ierr)
#endif

  return

end subroutine mp_send

subroutine mp_recv(val,n,tag,stat,ierr)

  implicit none

  real(8),intent(inout) :: val
  integer(4),intent(in) :: n,tag
  integer(4),intent(out) :: stat,ierr
   

#ifdef PARA
  call MPI_Recv(val,n,MPI_REAL,0,tag,MPI_COMM_WORLD,stat,ierr)
#endif

  return

end subroutine mp_recv

subroutine mp_reduce(val1,val2,n,ierr)

  implicit none

  real(8),intent(inout) :: val1,val2
  integer(4),intent(in) :: n
  integer(4),intent(out) :: ierr
  
#ifdef PARA
  call MPI_Reduce(val1,val2,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif

  return

end subroutine mp_reduce


subroutine mp_barrier(ierr)

  implicit none

  integer(4),intent(in) :: ierr

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

end subroutine mp_barrier

end module mp
