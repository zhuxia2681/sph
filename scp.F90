program scp
!
! Original code is icm.cpp
!            by G. Samsonidze.
!
! JDY: 2021.2.16
!

  use number
  use numerical
  use io
  use cell
  use integral
  use sli1
  use mp
  use util

  implicit none

!//JDY
  integer(4) :: i,j,k
  real(8),allocatable :: ravg(:,:,:)
  integer(4) :: ifile,nfile

!// read input
  character(7) :: tmpchar
  character(20) :: fname
  character(10) :: form
  logical :: uc,sc,renorm
  real(8) :: uco(3),sco(3),ucv(3,3),scv(3,3)
  real(8) :: thr
  integer(4) :: l1,l2,m1,m2
  integer(4) :: thrpow,coulpow
  integer(4) :: ucu,scu
  integer(4) :: rank,ierr,size
!JDY 2021.2.25 MPI
  integer(4) :: tag,stat,ii,jj,chkrank
  integer(4) :: flag1,flag2,flag3
!//
  real(8) :: vol,origin(3),span(3,3),stepvec(3,3)
  real(8) :: vol1,origin1(3),span1(3,3),stepvec1(3,3)
  real(8) :: vol2,origin2(3),span2(3,3),stepvec2(3,3)
  real(8) :: cellvec(3,3),cellvec1(3,3),cellvec2(3,3)
  integer(4) :: ngrid(3),ngrid1(3),ngrid2(3),tmpgrid(3)
  real(8),allocatable :: val(:,:,:),tmpval(:,:,:)
  real(8),allocatable :: val1(:,:,:)
  real(8),allocatable :: val2(:,:,:)

!//write_xsf
  character(20) :: tmpfname

!//
  real(8) :: tot

!// isovalue_scale
  real(8) :: isovalue(2)

!// scalar_trunc
  real(8) :: tmpori(3)
  real(8),allocatable :: tmpsf(:,:,:)

!// scalar_norm
  real(8),allocatable :: normedval(:,:,:)
  real(8) :: wfnorm

!// conv2bohr
  real(8) :: tmpvec(3,3), tmpo(3)

!// set_cutoff
  real(8) :: rcutoff,r2cutoff

!// check_overlap
  logical :: flag_overlap

!// grid_offset
  real(8) :: offset(3)

!// mc_average
  integer(4) :: nnn

!//
  real(8) :: tmpene,energy

!//ferdi
  integer(4) :: l,m
  complex(8),allocatable :: Rnl(:,:,:)

!// clock
  integer(4) :: t1,t2
  real(8) :: t

!//PARA
#ifdef PARA
  integer(4) :: idummy,info
#endif


#ifdef PARA
  ierr = 0
  call mp_ini(ierr)
  if (ierr /= MPI_SUCCESS) then
    write(*,*) " mpi ini. return error."
  end if
  call mp_info(rank,size,ierr)
#else
  rank = 0
  size = 1
#endif

  if (rank == 0) call clock_start(t1)


  allocate(ravg(2*NCELL+1,2*NCELL+1,2*NCELL+1))
  ravg = 0.0

  if (rank == 0) then
    open(2,file='scp.inp',form="formatted",status='old')
    open(8,file='scp.log',form="formatted",status='unknown')
    read(2,*) tmpchar,nfile

    do ifile = 1, nfile

      info = 0
      call par_read(ifile,fname,form,thr,thrpow,coulpow,uc,&
           ucu,sc,scu,renorm,info,l1,l2,m1,m2)
      if (info /= 0) then
        write(*,*) " error: par_read ifile =",ifile
      end if

      info = 0
      call par_read_cell(ifile,uco,ucv,sco,scv,info)
      if (info /= 0) then
        write(*,*) " error: par_read_cell ifile =",ifile
      end if

      info = 0
      call read_xsf(fname,ngrid,origin,span,val,stepvec,info)

      if (info /= 0) then
        write(*,*) " error: read_xsf ifile =",ifile
      end if

      allocate(Rnl(ngrid(1),ngrid(2),ngrid(3)))
      if (ifile == 1) then
        l = l1
        m = m1
      else
        l = l2
        m = m2
      end if
      call ferdi(l,m,val,ngrid,origin,stepvec,Rnl)
      val = real(Rnl)
      deallocate(Rnl)
      write(tmpfname,'(a4,i1,a4)') "rad_",ifile,".xsf"
      call write_xsf(tmpfname,origin,stepvec,ngrid,val)

      call trimming(ngrid,val,tmpgrid,tmpval)
      ngrid = tmpgrid
      val = tmpval

      info = 0
      call cell_set(ngrid,uc,ucu,origin,span,uco,ucv,cellvec)
      if (info /= 0) then
        write(*,*) " error: cell_set ifile =",ifile
      end if

      if (thr < (1.0 -EPS9)) then

        call isovalue_scale(thrpow,thr,isovalue(ifile),val,ngrid)
        call scalar_trunc(isovalue(ifile),ngrid,val,stepvec,&
                          origin,tmpsf,tmpori,info)
      else
        isovalue(ifile) = 0.0
      end if

      info = 0
      call scalar_norm(coulpow,renorm,ngrid,stepvec,tmpsf,&
                       normedval,wfnorm)

      call conv2bohr(tmpori,stepvec,origin,tmpvec)

      if (ifile == 1) then
        ngrid1 = ngrid
        origin1 = origin
        stepvec1 = tmpvec
        cellvec1 = cellvec
        val1 = normedval
      else
        ngrid2 = ngrid
        origin2 = origin
        stepvec2 = tmpvec
        cellvec2 = cellvec
        val2 = normedval
      end if
    end do
  end if
#ifdef PARA
!  if (rank == 0) then

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ngrid1,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ngrid2,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!debug
  write(*,*) " rank ",rank
  if (rank /= 0) then
    allocate(val1(ngrid1(1),ngrid1(2),ngrid1(3)))
    allocate(val2(ngrid2(1),ngrid2(2),ngrid2(3)))
  end if
!    val1 = 0.0
!    val2 = 0.0
    do i = 1, ngrid1(1)
      do j = 1, ngrid1(2)
        do k = 1, ngrid1(3)
          call MPI_Bcast(val1(i,j,k),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        end do
      end do
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_Bcast(l1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(l2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(m1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(m2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(origin1,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(origin2,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    do i = 1, 3
      call MPI_Bcast(stepvec1(i,:),3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    end do
    do i = 1, 3
      call MPI_Bcast(stepvec2(i,:),3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    end do
!debug
if (rank /= 0) then
  write(*,*) " stepvec1",stepvec1
  write(*,*) " stepvec2",stepvec2
end if
    do i = 1, ngrid2(1)
      do j = 1, ngrid2(2)
        do k = 1, ngrid2(3)
          call MPI_Bcast(val2(i,j,k),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        end do
      end do
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
!  end if
#endif

  if (rank == 0) then
      write(*,*)
      write(*,'(a23)') " Parameters for grid 1"
      write(*,'(a10,x23i4)') "  grid = ",ngrid1
      write(*,*) "  Origin and step vectors "
      write(*,'(a4,3f13.9,a8)') "  o1",origin1,"(Bohr)"
      do i = 1, 3
        write(*,'(a4,i1,a1,3f13.9,a8)') "  v[",i,"]",stepvec1(i,:),"(Bohr)"
      end do
      write(*,'(a15,1f13.6)') " wfn isovalue =",isovalue(1)
      write(*,'(a15,1f13.6)') " wfn norm     =",wfnorm
      write(*,*)
      write(*,'(a23)') " Parameters for grid 2"
      write(*,'(a10,x23i4)') "  grid = ",ngrid2
      write(*,*) "  Origin and step vectors "
      write(*,'(a4,3f13.9,a8)') "  o1",origin2,"(Bohr)"
      do i = 1, 3
        write(*,'(a4,i1,a1,3f13.9,a8)') "  v[",i,"]",stepvec2(i,:),"(Bohr)"
      end do
      write(*,'(a15,1f13.6)') " wfn isovalue =",isovalue(2)
      write(*,'(a15,1f13.6)') " wfn norm     =",wfnorm
      write(*,*)
  end if

  if (rank == 0) then
    if (nfile == 1) then
      call check_overlap(rcutoff,ngrid1,ngrid2,val1,val2,stepvec1,stepvec2,&
           origin1,origin2,flag_overlap)
    else
      flag_overlap = .true.
    end if
  end if

  if (rank == 0) then
    if (flag_overlap) then
      write(*,*) " Wavefunction overlaps with its image. &
 Using Monte Carlo averaging."
      call grid_offset(origin1,stepvec1,origin2,stepvec2,offset)
    end if
  end if


#ifdef PARA
  call MPI_Bcast(offset,3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif


  call mc_average(offset,stepvec1,stepvec2,ravg)

!JDY: 
!    each rank already have values of ravg?
!#ifdef PARA
!  call MPI_Bcast(ravg,(2*NCELL+1)**3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!#endif


  call coulomb_integral(rank,size,ngrid1,ngrid2,val1,val2,stepvec1,stepvec2,&
       origin1,origin2,ravg,energy,ierr)

  if (rank == 0) then
    write(*,'(a21,1f14.6,a3)') "  coulomb integral = ",energy," eV"
  end if


#ifdef PARA
  ierr = 0
  call MPI_Finalize(ierr)
  if (ierr /= 0) then
    write(*,*) " mpi fin. return error."
  end if
#endif

if (rank == 0) then
    call clock_end(t1,t)

    write(*,*) 
    write(*,'(a10,1f13.6,a4)') " end time ",t,"sec."
    write(*,*) 

endif

end program scp
