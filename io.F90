module io

  use number
  use mp

  implicit none

contains

subroutine par_read(ifile,fname,form,thr,thrpow,coulpow,uc,&
       ucu,sc,scu,renorm,info,l1,l2,m1,m2)

    implicit none
  
    integer(4),intent(in) :: ifile

    character(20),intent(out) :: fname
    character(10),intent(out) :: form
    logical,intent(out) :: uc,sc,renorm
    real(8),intent(out) :: thr
    integer(4),intent(out) :: thrpow,coulpow
    integer(4),intent(out) :: ucu,scu
    integer(4),intent(out) :: l1,l2,m1,m2,info
 
    character(30) :: tmp
    character(20) :: s1,s2
    integer(4) :: rank,size,ierr
  

  call mp_info(rank,size,ierr)

  rewind(2)

!JDY
!  if (rank == 0) then
    write(*,*) " read input parameters."
    write(*,*) " ifile",ifile
    if (ifile == 2) then
      do 
        read(2,*,end=999) s1,s2
        if (s1 == "l2") read(s2,*) l2
        if (s1 == "m2") read(s2,*) m2
        if (s1 == "inputfilename2") fname = s2
        if (s1 == "inputfileformat2") form = s2
        if (s1 == "threshold2") read(s2,*) thr
        if (s1 == "threshold_power2") read(s2,*) thrpow
        if (s1 == "coulomb_power2") read(s2,*) coulpow
        if (s1 == "uc2") then 
          call boolchk(s2,uc,ierr)
          if (ierr /= 0) then 
            write(*,*) " uc2 incorrect setting. "
            stop
          end if
        end if
        if (s1 == "ucu2") then 
          call unitchk(s2,ucu,ierr)
          if (ierr /= 0) then 
            write(*,*) " ucu2 incorrect setting. "
            stop
          end if
        end if
        if (s1 == "sc2") then 
          call boolchk(s2,sc,ierr)
          if (ierr /= 0) then 
            write(*,*) " sc2 incorrect setting. "
            stop
          end if
        end if
        if (s1 == "scu2") then 
          call unitchk(s2,scu,ierr)
          if (ierr /= 0) then 
            write(*,*) " scu2 incorrect setting. "
            stop
          end if
        end if
        if (s1 == "renormalize2") then 
          call boolchk(s2,renorm,ierr)
          if (ierr /= 0) then 
            write(*,*) " renorm2. incorrect setting. "
            stop
          end if
        end if
      end do
    elseif (ifile == 1) then
      do 
        read(2,*,end=999) s1,s2
        if (s1 == "l1") read(s2,*) l1
        if (s1 == "m1") read(s2,*) m1
        if (s1 == "inputfilename") fname = s2
        if (s1 == "inputfileformat") form = s2
        if (s1 == "threshold") read(s2,*) thr
        if (s1 == "threshold_power") read(s2,*) thrpow
        if (s1 == "coulomb_power") read(s2,*) coulpow
        if (s1 == "uc") then 
          call boolchk(s2,uc,ierr)
          if (ierr /= 0) then 
            write(*,*) " uc incorrect setting. "
            stop
          end if
        end if
        if (s1 == "ucu") then 
          call unitchk(s2,ucu,ierr)
          if (ierr /= 0) then 
            write(*,*) " ucu incorrect setting. "
            stop
          end if
        end if
        if (s1 == "sc") then 
          call boolchk(s2,sc,ierr)
          if (ierr /= 0) then 
            write(*,*) " sc incorrect setting. "
            stop
          end if
        end if
        if (s1 == "scu") then 
          call unitchk(s2,scu,ierr)
          if (ierr /= 0) then 
            write(*,*) " scu incorrect setting. "
            stop
          end if
        end if
        if (s1 == "renormalize") then 
          call boolchk(s2,renorm,ierr)
          if (ierr /= 0) then 
            write(*,*) " renorm. incorrect setting. "
            stop
          end if
        end if
      end do 
    end if
  999 write(*,*) " end read parameters."
!  end if
  info = 0

end subroutine par_read

subroutine par_read_cell(ifile,uco,ucv,sco,scv,info)

    implicit none
  
    integer(4),intent(in) :: ifile
    character(30) :: tmp
    real(8),intent(out) :: uco(3),sco(3),ucv(3,3),scv(3,3)
    integer(4),intent(out) :: info

    integer(4) :: i
    integer(4) :: rank,size,ierr
  

  call mp_info(rank,size,ierr)

!  if (rank == 0) then
    write(*,*) " read cell parameters."
    rewind(2)
    do 
      read(2,*,end=999) tmp
      if (ifile == 2) then
        if (trim(tmp) == "uco2") then
          read(2,*) uco(1),uco(2),uco(3)
        end if
        if (trim(tmp) == "ucv2") then
          do i=1,3
            read(2,*) ucv(1,i),ucv(2,i),ucv(3,i)
          end do
        end if
        if (trim(tmp) == "sco2") then
          read(2,*) sco(1),sco(2),sco(3)
        end if
        if (trim(tmp) == "scv2") then
          do i=1,3
            read(2,*) scv(1,i),scv(2,i),scv(3,i)
          end do
        end if
      else
        if (trim(tmp) == "uco") then
          read(2,*) uco(1),uco(2),uco(3)
        end if
        if (trim(tmp) == "ucv") then
          do i=1,3
            read(2,*) ucv(1,i),ucv(2,i),ucv(3,i)
          end do
        end if
        if (trim(tmp) == "sco") then
          read(2,*) sco(1),sco(2),sco(3)
        end if
        if (trim(tmp) == "scv") then
          do i=1,3
            read(2,*) scv(1,i),scv(2,i),scv(3,i)
          end do
        end if
      end if
    end do

999 write(*,*) " end cell parameters."
!  end if

  info = 0

end subroutine par_read_cell


subroutine boolchk(boolin,boolout,ierr)

  implicit none

  character(10),intent(in) :: boolin
  logical,intent(out) :: boolout
  integer(4),intent(out) :: ierr

  if (trim(boolin) == "T" .or. trim(boolin) == "True" .or. &
      trim(boolin) == "true" .or. trim(boolin) == "t" .or. &
      trim(boolin) == ".true.") then
    boolout = .true.
    ierr = 0
  elseif (trim(boolin) == "F" .or. trim(boolin) == "False" .or. &
      trim(boolin) == "false" .or. trim(boolin) == "f" .or. &
      trim(boolin) == ".false.") then
    boolout = .false.
    ierr = 0
  else
    ierr = -1
  end if

end subroutine boolchk

subroutine unitchk(unittypein,unittypeout,ierr)

  implicit none

  character(10),intent(in) :: unittypein
  integer(4),intent(out) :: unittypeout,ierr

  if (trim(unittypein) == "bohr") then
    unittypeout = 0
    ierr = 0
  elseif (trim(unittypein) == "angstrom") then
    unittypeout = 1
    ierr = 0
  elseif (trim(unittypein) == "latvec") then
    unittypeout = 2
    ierr = 0
  else
    ierr = -1
  end if

end subroutine unitchk

subroutine read_xsf(fname,ngrid,origin,span,val,stepvec,info)
!
! JDY: 2021.2.15
!

  implicit none

  character(20),intent(in) :: fname
  integer(4),intent(out) :: ngrid(3)
  real(8),intent(out) :: origin(3),span(3,3),stepvec(3,3)
  real(8),allocatable :: tmpval(:,:,:)
  real(8),allocatable,intent(out) :: val(:,:,:)
  integer(4),intent(out) :: info
!//
  character(30) :: tmp
  integer(4) :: i,j,k,ix,iy,iz,tmpgrid(3)
  real(8) :: tmpvec(3,3)
  integer(4) :: rank,size,ierr,flag

!//debug
  real(8) :: debug

  rewind(11)

  call mp_info(rank,size,ierr)
  
  open(11,file=fname,form='formatted',status='old')

!debug
  write(*,*) " read in "
 
  do 
    read(11,*,end=999) tmp
    if (trim(tmp) == "BEGIN_BLOCK_DATAGRID_3D") then
      read(11,*,end=999) 
      read(11,*,end=999) 
      read(11,*) (tmpgrid(i),i=1,3)
      read(11,*) (origin(i),i=1,3)
      do j = 1, 3
        read(11,*) (span(i,j),i=1,3)
      end do

      allocate(tmpval(tmpgrid(1),tmpgrid(2),tmpgrid(3)))
      read(11,*) (((tmpval(ix,iy,iz),ix=1,tmpgrid(1)),iy=1,tmpgrid(2)),&
                  iz=1,tmpgrid(3))
      read(11,*,end=999) tmp
      if (tmp /= "END_DATAGRID_3D") then
        ierr =-1
      else
        ierr = 0
      end if
    end if
  end do 

999 write(*,*) " read xsf ",trim(fname)," end."

!JDY: 2021.2.22
!JDY: 2021.2.24 comment out.
  do i = 1, 3
    ngrid(i) = tmpgrid(i) !- 1
  end do
!end

  allocate(val(ngrid(1),ngrid(2),ngrid(3)))
  do ix = 1, ngrid(1)
    do iy = 1, ngrid(2)
      do iz = 1, ngrid(3)
        val(ix,iy,iz) = tmpval(ix,iy,iz)
      end do
    end do
  end do

!JDY: compute step vectors.
  do i = 1, 3
    tmpvec(i,i) = origin(i)
  end do

  span  = span - tmpvec

  do i = 1, 3
    do j = 1, 3
      stepvec(i,j) = span(i,j) / dble(ngrid(j))
    end do
  end do

  close(11)
  info = 0

end subroutine read_xsf

subroutine write_xsf(fname,origin,stepvec,ngrid,val)

  implicit none

  character(20),intent(in) :: fname
  character(20) :: tmpname
  real(8),intent(in) :: origin(3),stepvec(3,3)
  integer(4),intent(in) :: ngrid(3)
  real(8),allocatable :: val(:,:,:)
!local
  integer(4) :: i,j,k
  
  tmpname = trim(fname)
  open(100,file=tmpname,form='formatted',status='replace')

  write(100,*) "#"
  write(100,*) "# Generated by SCP."
  write(100,*) "#"
  write(100,*) "BEGIN_BLOCK_DATAGRID_3D"
  write(100,*) "3D_field"
  write(100,*) "BEGIN_DATAGRID_3D_UNKNOWN"
  write(100,'(3i6)') ngrid
  write(100,'(3f13.6)') origin
  write(100,'(3f13.9,/,3f13.9,/,3f13.9)') stepvec
  write(100,'(6f13.9)') (((val(i,j,k),i=1,ngrid(1)),j=1,ngrid(2)),&
   k=1,ngrid(3))  
  write(100,*) "END DATAGRID_3D"
  write(100,*) "END_BLOCK_DATAGRID_3D"
 
  close(100)

end subroutine write_xsf

end module io
