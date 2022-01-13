module integral

  use number
  use numerical
  use cell
  use mp

  implicit none

  contains

subroutine scalar_trunc(thr,ngrid,sf,stepvec,sfoin,sfout,sfoout,info)

  implicit none

  real(8),intent(in) :: thr,sfoin(3),stepvec(3,3)
  real(8),allocatable,intent(in) :: sf(:,:,:)
  real(8),allocatable,intent(out) :: sfout(:,:,:)
  real(8),intent(out) :: sfoout(3)
  integer(4),intent(in) :: ngrid(3)
  integer(4),intent(out) :: info

!JDY local
  integer(4) :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  integer(4) :: ni,nj,nk,ni1,nj1,nk1,ni2,nj2,nk2
  real(8),allocatable :: sf0(:,:,:)

  ni1 = ngrid(1)
  nj1 = ngrid(2)
  nk1 = ngrid(3)
  
  imin = ni1
  imax = -1
  jmin = nj1
  jmax = -1
  kmin = nk1
  kmax = -1

  do i = 1, ni1
    do j = 1, nj1
      do k = 1, nk1
        if (abs(sf(i,j,k)) > thr) then
          if (i < imin) imin = i 
          if (i > imax) imax = i 
          if (j < jmin) jmin = j 
          if (j > jmax) jmax = j 
          if (k < kmin) kmin = k 
          if (k > kmax) kmax = k 
        end if
      end do
    end do
  end do

  ni2 = imax - imin + 1
  nj2 = jmax - jmin + 1
  nk2 = kmax - kmin + 1

  if (ni2 < 1 .or. nj2 < 1 .or. nk2 < 1) then
    write(*,*) " Error in scalar_trunc",ni2,nj2,nk2
    stop
  end if


  if (ni2 /= ni1 .or. nj2 /= nj1 .or. nk2 /= nk1) then
    ni = ni2
    nj = nj2
    nk = nk2
    sf0 = sf

    do i = 1, ni2 
      do j = 1, nj2 
        do k = 1, nk2 
          sfout(i,j,k) = sf0(i-1+imin, j-1+jmin, k-1+kmin) 
        end do
      end do
    end do

    sfoout(1) = sfoin(1) + stepvec(1,1)*dble(imin-1) + stepvec(2,1)*dble(imin-1)+ stepvec(3,1)*dble(kmin-1)
    sfoout(2) = sfoin(2) + stepvec(1,2)*dble(imin-1) + stepvec(2,2)*dble(jmin-1)+ stepvec(3,2)*dble(kmin-1)
    sfoout(3) = sfoin(3) + stepvec(1,3)*dble(imin-1) + stepvec(2,3)*dble(jmin-1)+ stepvec(3,3)*dble(kmin-1)
  end if

  if (ni2 == ni1 .or. nj2 == nj1 .or. nk2 == nk1) then
    write(*,*) " WARNING: The selected isovalue overlaps with&
 the edge of the cell. "
    write(*,*) " Beware, the Coulomb integral does not use periodic &
 boundary conditions."
  sfout = sf
  sfoout = sfoin
  else
    info = 0
  end if

end subroutine scalar_trunc

subroutine scalar_norm(coulomb_power,renorm,ngrid,stepvec,sf,&
           sfout,wfnorm)

  implicit none

  integer(4),intent(in) :: coulomb_power,ngrid(3)
  real(8),allocatable,intent(in) :: sf(:,:,:)
  real(8),allocatable,intent(out) :: sfout(:,:,:)
  logical,intent(in) :: renorm
  real(8),intent(in) :: stepvec(3,3)
  real(8),intent(out) :: wfnorm

!JDY local
  integer(4) :: i,j,k,ierr
  real(8) :: vv,ww,vol,ww2
  real(8),allocatable :: sftmp(:,:,:)

  ierr = 0
  ww = 0.0
  allocate(sftmp(ngrid(1),ngrid(2),ngrid(3)))
  allocate(sfout(ngrid(1),ngrid(2),ngrid(3)))
  sftmp = 0.0
  sfout = 0.0

  do i = 1, ngrid(1)
    do j = 1, ngrid(2)
      do k = 1, ngrid(3)
        if (coulomb_power == 1) then
          sftmp(i,j,k) = abs(sf(i,j,k)) 
        else
          sftmp(i,j,k) = sf(i,j,k)**2
        end if
      end do
    end do
  end do

  call volume(stepvec,vol)

  vv = vol / BOHR**3

  if (vv < EPS9) ierr = -1

  do i = 1, ngrid(1)
    do j = 1, ngrid(2)
      do k = 1, ngrid(3)
        ww = ww + sftmp(i,j,k)
      end do
    end do
  end do

  if (ww < EPS9) ierr = -1

  ww = ww * vv

  wfnorm = ww

  if (renorm) then
    do i = 1, ngrid(1) 
      do j = 1, ngrid(2) 
        do k = 1, ngrid(3) 
          sftmp(i,j,k) = sftmp(i,j,k) / wfnorm
        end do
      end do
    end do
  end if

  sfout = sftmp
  deallocate(sftmp)

end subroutine scalar_norm


subroutine check_overlap(rcutoff,FFTgrid1,FFTgrid2,sf1,sf2,stepvec1,stepvec2,&
           sf1o,sf2o,flag_overlap)

  implicit none

  real(8),intent(in) :: rcutoff
  integer(4) :: rank,size,ierr
  integer(4) :: FFTgrid1(3),FFTgrid2(3)
  real(8),allocatable,intent(in) :: sf1(:,:,:),sf2(:,:,:)
  real(8) :: stepvec1(3,3),stepvec2(3,3),sf1o(3),sf2o(3)

  integer(4) :: i,nplus,nminus,iminus,ngrid1,ngrid2
  integer(4) :: dplus,dminus
  real(8) :: vl1(3),vl2(3),vdotv,vcosv,vdotr,pp(4),dd(4)
  real(8) :: r(4,3)

  logical,intent(out) :: flag_overlap

#ifdef PARA 
  call mp_info(rank,size,ierr) 
#else
  rank = 0
  size = 1
#endif

  do i = 1, 3
    vl1(i) = sqrt(stepvec1(i,1)**2 + stepvec1(i,2)**2 + stepvec1(i,3)**2)
    vl2(i) = sqrt(stepvec2(i,1)**2 + stepvec2(i,2)**2 + stepvec2(i,3)**2)
    vdotv = stepvec1(i,1) * stepvec2(i,1) + stepvec1(i,2) * stepvec2(i,2) + stepvec1(i,3) * stepvec2(i,3)
    vcosv = vdotv / (vl1(i) * vl2(i))            

    if (abs(vcosv - 1.0) < EPS9) nplus = nplus + 1
    if (abs(vcosv + 1.0) < EPS9) then
      nminus = nminus + 1
      iminus = i
    end if
  end do

  if (nplus == 2 .and. nminus == 1) then
    if (iminus == 0) then
      ngrid1 = FFTgrid1(1)
      ngrid2 = FFTgrid2(1)
    elseif (iminus == 1) then
      ngrid1 = FFTgrid1(2)
      ngrid2 = FFTgrid2(2)
    elseif (iminus == 2) then
      ngrid1 = FFTgrid1(3)
      ngrid2 = FFTgrid2(3)
    end if
 
    r = 0.0
    r(1,1) = sf1o(1) 
    r(1,2) = sf1o(1) 
    r(1,3) = sf1o(2) 
    r(2,1) = sf1o(1) + stepvec1(iminus,1) * dble(ngrid1 - 1) 
    r(2,2) = sf1o(2) + stepvec1(iminus,2) * dble(ngrid1 - 1) 
    r(2,3) = sf1o(3) + stepvec1(iminus,3) * dble(ngrid1 - 1) 

    r(3,1) = sf2o(1) 
    r(3,2) = sf2o(1) 
    r(3,3) = sf2o(2) 
    r(4,1) = sf2o(1) + stepvec2(iminus,1) * dble(ngrid2 - 1) 
    r(4,2) = sf2o(2) + stepvec2(iminus,2) * dble(ngrid2 - 1) 
    r(4,3) = sf2o(3) + stepvec2(iminus,3) * dble(ngrid2 - 1) 

    do i = 1, 4
      vdotr = stepvec1(iminus,1) * r(i,1) + stepvec1(iminus,2) * r(i,2) + stepvec1(iminus,3) * r(i,3) 
      pp(i) = vdotr / vl1(iminus)
    end do 

    dd(1) = pp(1) - pp(3)
    dd(2) = pp(1) - pp(4)
    dd(3) = pp(2) - pp(3)
    dd(4) = pp(2) - pp(4)

    do i = 1, 4
      if (dd(i) > 0.0) then
        dplus = dplus + 1
      else
        dminus = dminus + 1
      end if
    end do

    if (dplus ==4 .or. dminus == 4) then
      do i = 1, 4
        if (abs(dd(i)) < rcutoff) then
          flag_overlap = .true.
        end if
      end do
    else
      flag_overlap = .true.
    end if

  else
    if (rank == 0) then
      write(*,*) " The image plane is not parallel/perpendicular&
 to the lattice vectors."
      write(*,*) " Skipping wavefunction overlap check and Monte&
 Carlo averaging."
      write(*,*) " Your job may fail if the wavefunction overlaps &
  with its image."
    end if    
  end if

end subroutine check_overlap


subroutine grid_offset(sfo1,stepvec1,sfo2,stepvec2,offset)

  implicit none

  real(8),intent(in) :: sfo1(3),sfo2(3)
  real(8),intent(in) :: stepvec1(3,3),stepvec2(3,3)
  integer(4) :: i,j,nn(3)
  real(8) :: lv2,vdotd,dd(3)
  real(8),intent(out) :: offset(3)

  do i = 1,3
    dd(i) = sfo1(i) - sfo2(i) 
  end do

  do i = 1, 3
    lv2 = stepvec1(i,1)**2 + stepvec1(i,2)**2 + stepvec1(i,3)**2
    vdotd = stepvec1(i,1)*dd(1) + stepvec1(i,2)*dd(2) + stepvec1(i,3)*dd(3)

    call double_to_int(vdotd/lv2, nn(i))

    dd(1) = dd(1) - stepvec1(i,1)*dble(nn(1))
    dd(2) = dd(2) - stepvec1(i,2)*dble(nn(2))
    dd(3) = dd(3) - stepvec1(i,3)*dble(nn(3))
  end do

  offset = dd

end subroutine grid_offset


subroutine mc_average(offset,v1,v2,ravg)

  implicit none

  integer(4) :: rank,size,ierr
  real(8),intent(in) :: offset(3),v1(3,3),v2(3,3)
  real(8),allocatable,intent(out) :: ravg(:,:,:)
! local
  integer(4) :: i,j,k,l,m,lmin,lmax,navg
  integer(4) :: nlocal(2*NCELL+1,2*NCELL+1,2*NCELL+1)
  
#ifdef PARA
  integer(4) :: ndummy(2*NCELL+1,2*NCELL+1,2*NCELL+1)
#endif 

  real(8) :: rinv,r2,rr(3)
  real(8) :: rlocal(2*NCELL+1,2*NCELL+1,2*NCELL+1)

#ifdef PARA
  real(8) :: rdummy(2*NCELL+1,2*NCELL+1,2*NCELL+1)
#endif 
  
  real(8) :: dd(3),crand(3)
!JDY
  real(8) :: tmp
  integer(4) :: pow

#ifdef PARA
  call mp_info(rank,size,ierr)
  pow = (2*NCELL+1)**3
#else
  rank = 0
  size = 1
#endif

  allocate(ravg(2*NCELL+1,2*NCELL+1,2*NCELL+1))
  ravg = 0.0

#ifdef VERBOSE
  if (rank == 0) then
    write(*,*) " averaging 1 / | r - r' | over r' in the grid cell."
  end if
#endif

  call double_to_int(dble(NRAND) * dble(rank) / dble(size),lmin) 
  call double_to_int(dble(NRAND) * dble(rank+1) / dble(size),lmax) 

  do i = -NCELL,NCELL
    do j = -NCELL,NCELL
      do k = -NCELL,NCELL
        nlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) = lmax - lmin
        rlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) = 0.0
      end do
    end do
  end do

  do l = lmin, lmax
    do m = 1, 3
      call rand(tmp)  
      rr(m) = tmp - 0.5
    end do

    crand(1) = v2(1,1) * rr(1) + v2(2,1) * rr(2) + v2(3,1) * rr(3)
    crand(2) = v2(1,2) * rr(1) + v2(2,2) * rr(2) + v2(3,2) * rr(3)
    crand(3) = v2(1,3) * rr(1) + v2(2,3) * rr(2) + v2(3,3) * rr(3)

    do i = -NCELL,NCELL
      do j = -NCELL,NCELL
        do k = -NCELL,NCELL
          dd(1) = offset(1) + v1(1,1) + dble(i) + v1(2,1) + dble(j) + v1(3,1) + dble(k)
          dd(2) = offset(2) + v1(1,2) + dble(i) + v1(2,2) + dble(j) + v1(3,2) + dble(k)
          dd(3) = offset(3) + v1(1,3) + dble(i) + v1(2,3) + dble(j) + v1(3,3) + dble(k)
          r2 = (dd(1) - crand(1))**2 + (dd(2) - crand(2))**2 + (dd(3) - crand(3))**2 
          if (r2 > EPS9) then
            rlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) = rlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) + 1.0/sqrt(r2) 
          else
            nlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) = nlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) - 1 
          end if
        end do
      end do
    end do
  end do

#ifdef PARA
  do i = 1, size
    call MPI_ALLREDUCE(nlocal,ndummy,pow,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

    if (ierr /= 0) then
      write(*,*) " error: allreduce failed."
      stop
    end if
  end do

  do i = -NCELL,NCELL
    do j = -NCELL,NCELL
      do k = -NCELL,NCELL
        nlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) = ndummy(i+NCELL+1,j+NCELL+1,k+NCELL+1)
      end do
    end do
  end do

  call MPI_ALLREDUCE(rlocal,rdummy,pow,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  
    
  do i = -NCELL,NCELL
    do j = -NCELL,NCELL
      do k = -NCELL,NCELL
        rlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1) = rdummy(i+NCELL+1,j+NCELL+1,k+NCELL+1)
      end do
    end do
  end do
#endif

  do i = -NCELL,NCELL
    do j = -NCELL,NCELL
      do k = -NCELL,NCELL
        ravg(i+NCELL+1,j+NCELL+1,k+NCELL+1) = rlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1)/dble(nlocal(i+NCELL+1,j+NCELL+1,k+NCELL+1))
      end do
    end do
  end do

end subroutine mc_average


subroutine tab_average(v1,p1,p2,offset,ravg,rinv)

  implicit none

  real(8),intent(in) :: p1(3),p2(3),v1(3,3),offset(3)

  logical :: f_offset,f_range
  integer(4) :: i,j,ii, nn(3)
  real(8) :: lv2,vdotd
  real(8) :: dd(3)
  real(8),intent(in) :: ravg(:,:,:)
  real(8),intent(out) :: rinv
!debug
  integer(4) :: rank,size,ierr

#ifdef PARA
  call mp_info(rank,size,ierr)
#endif

  dd = 0.0
  do i = 1, 3 
    dd(i) = p1(i) - p2(i) - offset(i) 
  end do 

  
  nn = 0
  do i = 1, 3
    lv2 = v1(i,1)**2 + v1(i,2)**2 + v1(i,3)**2
    vdotd = v1(i,1)*dd(1) + v1(i,2)*dd(2) + v1(i,3)*dd(3)
    call double_to_int((vdotd/lv2),nn(i))
   
    dd(1) = dd(1) - v1(i,1)*dble(nn(1))
    dd(2) = dd(2) - v1(i,2)*dble(nn(2))
    dd(3) = dd(3) - v1(i,3)*dble(nn(3))
  end do


  if (abs(dd(1)) .lt. EPS9 .and. abs(dd(2)) .lt. EPS9 .and. abs(dd(3)) .lt. EPS9) then
    f_offset = .true.
  end if
  if (abs(nn(1)) .le. NCELL .and. abs(nn(2)) .le. NCELL .and. abs(nn(3)) .le. NCELL) then
    f_range = .true.
  end if


  if (f_offset .and. f_range) then
    rinv = ravg(nn(1)+NCELL,nn(2)+NCELL,nn(3)+NCELL)
  elseif (.not. f_offset) then
    write(*,*) "Error: f_offset failed in tab_average."
    rinv = -1.0
  elseif (.not. f_range) then
    write(*,*) "Error: f_range failed in tab_average."
    rinv = -1.0
  end if


end subroutine tab_average


subroutine coulomb_integral(rank,size,ngrid1,ngrid2,sf1,sf2,stepvec1,stepvec2,&
            origin1,origin2,ravg,energy,ierr)

  implicit none

  integer(4),intent(in) :: rank,size
  real(8),allocatable,intent(in) :: sf1(:,:,:),sf2(:,:,:)
  real(8),intent(in) :: stepvec1(3,3),stepvec2(3,3),origin1(3),origin2(3)
  integer(4),intent(in) :: ngrid1(3),ngrid2(3)

  integer(4) :: i,j,k
  integer(4) :: i1,j1,k1,i2,j2,k2,nx,ny,nz
  real(8),allocatable :: pnl1(:,:,:),pnl2(:,:,:)

  real(8) :: vol1,vol2,p1(3),p2(3),r2,rinv,ww
  real(8) :: rcutoff,r2cutoff,offset(3)
  real(8),allocatable,intent(in) :: ravg(:,:,:)
  real(8),intent(out) :: energy

  integer(4) :: info

#ifdef PARA
  real(8) :: rdummy
#endif

#ifdef VERBOSE
  integer(4) :: pnew, pold
#endif

  integer(4),intent(out) :: ierr

!debug
  integer(4) :: cnt


!JDY step-volume check
  call volume(stepvec1,vol1)
  call volume(stepvec2,vol2)

  if (vol1 < EPS9) then
    write(*,*) " Error: cell volume 1 = 0 in coulomb_integral. "
    stop
  elseif (vol2 < EPS9) then
    write(*,*) " Error: cell volume 2 = 0 in coulomb_integral. "
    stop
  end if


  p1 = 0.0
  p2 = 0.0
#ifdef VERBOSE
  pnew = 0
  pold = 0
#endif

!JDY: 2021.2.21
!     how about ngrid2.
  nx = ngrid1(1)   
  ny = ngrid1(2)  
  nz = ngrid1(3)  

  do i1 = 1, nx 
    do j1 = 1, ny 
      do k1 = 1, nz 
#ifdef PARA
        if (mod((i1-1)*ny*nz + (j1-1)*nz + (k1-1),size) /= rank) then 
          cycle
        end if
#endif

#ifdef VERBOSE
        if (rank == 0) then
          call double_to_int((100.0*dble((i1-1)*ny*nz + (j1-1)*nz + (k1-1))) / dble(nx*ny*nz),pnew)
          if (pnew /= pold) then
            write(*,*) " completed ",pnew," %"
            pold = pnew
          endif
        endif
#endif
        p1(1) = origin1(1) + stepvec1(1,1)*dble(i1) + stepvec1(1,2)*dble(j1) + stepvec1(1,3)*dble(k1)
        p1(2) = origin1(2) + stepvec1(2,1)*dble(i1) + stepvec1(2,2)*dble(j1) + stepvec1(2,3)*dble(k1)
        p1(3) = origin1(3) + stepvec1(3,1)*dble(i1) + stepvec1(3,2)*dble(j1) + stepvec1(3,3)*dble(k1)

        do i2 = 1, nx
          do j2 = 1, ny
            do k2 = 1, nz
              p2(1) = origin2(1) + stepvec2(1,1)*dble(i2) + stepvec2(1,2)*dble(j2) + stepvec2(1,3)*dble(k2)
              p2(2) = origin2(2) + stepvec2(2,1)*dble(i2) + stepvec2(2,2)*dble(j2) + stepvec2(2,3)*dble(k2)
              p2(3) = origin2(3) + stepvec2(3,1)*dble(i2) + stepvec2(3,2)*dble(j2) + stepvec2(3,3)*dble(k2)
              r2 = (p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2
              if (r2 > r2cutoff) then
                rinv = 1.0 / sqrt(r2)   
              else
                call tab_average(stepvec1,p1,p2,offset,ravg,rinv)
              end if
              if (rinv < -EPS9) then
                write(*,*) "Error: rinv < 0 in coulomb_integral."
                ierr = -1
                stop
              end if

              ww = ww + sf1(i1,j1,k1) * sf2(i2,j2,k2) * rinv
!debug
              cnt = cnt + 1

            end do  
          end do  
        end do  
      end do
    end do
  end do

!debug
write(8,*) " cnt",cnt

#ifdef PARA
  call mp_barrier(ierr)
#endif


  ww = ww * 0.5 * HARTREE * vol1 * vol2


#ifdef PARA
  call mp_reduce(ww,rdummy,1,ierr)
  ww = rdummy
#endif

  energy = ww

end subroutine coulomb_integral


subroutine set_cutoff(v1,v2,rcutoff,r2cutoff)

  implicit none

  real(8),intent(in) :: v1(3,3),v2(3,3)
  real(8),intent(out) :: rcutoff,r2cutoff
  real(8) :: l2,l2min
  integer(4) :: ii

  l2min = INF9
  
  do ii = 1,3
    l2 = v1(ii,1)**2 + v1(ii,2)**2 + v1(ii,3)**2
    if (l2min > l2) l2min = l2
    l2 = v2(ii,1)**2 + v2(ii,2)**2 + v2(ii,3)**2
    if (l2min > l2) l2min = l2
  end do

  r2cutoff = l2min * dble(NCELL**2)
  rcutoff = sqrt(r2cutoff)

#ifdef PARA
  call MPI_Bcast(r2cutoff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD)
  call MPI_Bcast(rcutoff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD)
#endif

end subroutine set_cutoff

subroutine total(ngrid,val,tot)

  implicit none

  integer(4),intent(in) :: ngrid(3)
  real(8),allocatable,intent(in) :: val(:,:,:)
  real(8) :: ww
  real(8),intent(out) :: tot


  integer(4) :: i,j,k,nx,ny,nz

  nx = ngrid(1)
  ny = ngrid(2)
  nz = ngrid(3)
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        ww = ww + abs(val(i,j,k))
      end do
    end do
  end do

  tot = ww

end subroutine total


end module integral
