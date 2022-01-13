module symmetries

  use number

  implicit none

  contains

subroutine find_group(nrot,smat,gname,code_group)
!
!  Given a group of nrot rotation matrices omat (in cartesian coordinates)
!  this routine finds the name of the point group. It assumes but does not
!  check that:
!  1) The nrot matrices smat are actually a group.
!  2) The group is one of the thirty-two point groups.
!
  implicit none

  integer(4),intent(in) :: nrot
  integer(4),intent(out) :: code_group
  real(8),intent(in) :: smat(3,3,nrot)
  character(11) :: gname
  integer(4) :: noperation(6), irot, ts  
!
! For each possible group operation the function tipo_sym gives a code
!   1 identity,
!   2 inversion,
!   3 proper rotation <> 180,
!   4 proper rotation 180 degrees,
!   5 mirror,
!   6 improper rotation
! the variable noperation counts how many operations are present in the group.
!
  noperation = 0
  
  do irot=1,nrot
    ts = tipo_sym(smat(:,:,irot))
    noperation(ts) = noperation(ts) + 1
  end do

  if (noperation(1).ne.1) then
    write(*,*) "find_group the group has not identity E."
    stop
  endif

  code_group = 0
!
! table of code group

!  1 C_1
!  2 C_i
!  3 C_s
!  4 C_2
!  5 C_3
!  6 C_4
!  7 C_6
!  8 D_2
!  9 D_3
! 10 D_4
! 11 D_6
! 12 C_2v
! 13 C_3v
! 14 C_4v
! 15 C_6v
! 16 C_2h
! 17 C_3h
! 18 C_4h
! 19 C_6h
! 20 D_2h
! 21 D_3h
! 22 D_4h
! 23 D_6h
! 24 D_2d
! 25 D_3d
! 26 S_4
! 27 S_6
! 28 T
! 29 T_h
! 30 T_d
! 31 O
! 32 O_h
  
  if (noperation(2) == 1) then
!
!  There is not inversion
!
   IF (nrot==1) THEN
      code_group=1                                          ! C_1
   ELSEIF (nrot==2) THEN
      IF (noperation(4)==1) code_group=4                    ! C_2
      IF (noperation(5)==1) code_group=3                    ! C_s
   ELSEIF (nrot==3) THEN
      IF (noperation(3)==2) code_group=5                    ! C_3
   ELSEIF (nrot==4) THEN
      IF (noperation(6)>0)  code_group=26                   ! S_4
      IF (noperation(5)>0.and.code_group==0) code_group=12  ! C_2v
      IF (noperation(3)>0.and.code_group==0) code_group=6   ! C_4
      IF (noperation(4)>0.and.code_group==0) code_group=8   ! D_2
   ELSEIF (nrot==6) THEN
      IF (noperation(5)==3) code_group=13                   ! C_3v
      IF (noperation(5)==1) code_group=17                   ! C_3h
      IF (noperation(4)==3.and.code_group==0) code_group=9  ! D_3
      IF (noperation(3)>0.and.code_group==0) code_group=7   ! C_6
   ELSEIF (nrot==8) THEN
      IF (noperation(5)==4) code_group=14                   ! C_4v
      IF (noperation(5)==2) code_group=24                   ! D_2d
      IF (noperation(3)>0.and.code_group==0) code_group=10  ! D_4
   ELSEIF (nrot==12) THEN
      IF (noperation(5)==6) code_group=15                   ! C_6v
      IF (noperation(5)==4) code_group=21                   ! D_3h
      IF (noperation(4)>6.and.code_group==0) code_group=11  ! D_6
      IF (noperation(3)>0.and.code_group==0) code_group=28  ! T
   ELSEIF (nrot==24) THEN
      IF (noperation(5)>0) code_group=30                    ! T_d
      IF (noperation(5)==0) code_group=31                   ! O
   ELSE
      write(*,*) "find_group wrong number of elements"
      stop
   ENDIF
  ELSEIF (noperation(2)==1) THEN
  
  !
!  There is inversion
!
   IF (nrot==2) THEN
      code_group=2                                          ! C_i
   ELSEIF (nrot==4) THEN
      code_group=16                                         ! C_2h
   ELSEIF (nrot==6) THEN
      code_group=27                                         ! S_6
   ELSEIF (nrot==8) THEN
      IF (noperation(5)==3) code_group=20                   ! D_2h
      IF (noperation(5)==1) code_group=18                   ! C_4h
   ELSEIF (nrot==12) THEN
      IF (noperation(5)==3) code_group=25                   ! D_3d
      IF (noperation(5)==1) code_group=19                   ! C_6h
   ELSEIF (nrot==16) THEN
      IF (noperation(5)==5) code_group=22                   ! D_4h
   ELSEIF (nrot==24) THEN
      IF (noperation(5)>6) code_group=23                    ! D_6h
      IF (noperation(5)==3) code_group=29                   ! T_h
   ELSEIF (nrot==48) THEN
      code_group=32                                         ! O_h
   ELSE
    write(*,*) "find_group wrong number of elements"
    stop
   ENDIF
  ELSE
    write(*,*) "find_group too many inversions"
    stop
  ENDIF

  IF (code_group==0) then
    write(*,*) "'find_group incompatible operations"
    stop
  end if

  gname=group_name(code_group)

  RETURN

end subroutine find_group

subroutine divide_class(code_group,nrot,smat,nclass,nelem,elem,which_irr)
!
! This subroutine receives as input a set of nrot 3x3 matrices smat, which
! are assumed to be the operations of the point group given by code_group.
! smat are in cartesian coordinates.
! This routine divides the group in classes and find:
!
! nclass         the number of classes of the group
! nelem(iclass)  for each class, the number of elements of the class
! elem(i,iclass) 1<i<nelem(iclass) for each class tells which matrices
!                smat belong to that class
! which_irr(iclass) for each class gives the position of that class in the
!                character table associated with the group and provided
!                by the routine set_irr_rap. NB: changing the order of
!                the elements in the character table must be reflected in
!                a change to which_irr. Presently the character tables
!                are those given by P.W. Atkins, M.S. Child, and
!                C.S.G. Phillips, "Tables for group theory".
!                Several equivalent names for the irreducible representation
!                are given. D, G, L, S are used for Delta, Gamma, Lambda
!                and Sigma.
!
  implicit none

  integer(4)  :: &
          code_group,  &   ! The code of the point group
          nrot,        &   ! The number of symmetry operation
          nclass,      &   ! The number of classes
          nelem(12),   &   ! The elements of each class
          elem(8,12),  &   ! Which elements in the smat list for each class
          which_irr(12)    ! See above

  real(8) :: smat(3,3,nrot), cmat(3,3), ax(3), bx(3), cx(3), ars
  integer(4) :: ipol, axis, axis1, axis2, ts, iax, ibx, icx, aclass, &
           bclass, cclass, imax, imbx, imcx, amclass, bmclass, cmclass, ind2(3)
 
  real(8) :: angle_vectors, ax_save(3,2:4)

  logical :: is_parallel, isok, isok1, done_ax(6)

  eps = EPS7

!
! Divide the group in classes.
!
  nclass=0
  nelem=0
  done=0
  DO irot=1,nrot
    IF (done(irot)==0) THEN
      nclass=nclass+1
      DO jrot=1,nrot
        CALL coniug_mat(smat(1,1,jrot),smat(1,1,irot),cmat)
        DO krot=1,nrot
          IF (compare_mat(cmat,smat(1,1,krot)).AND.done(krot)==0) THEN
            nelem(nclass)=nelem(nclass)+1
            elem(nelem(nclass),nclass)=krot
            done(krot)=1
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDDO

!
!  For each class we should now decide which_irr. This depends on the group
!  and on the tables of characters of the irreducible representations,
!  so we must make different things for different groups.
!

  which_irr(1)=1
  IF (code_group==1) THEN
    IF (nclass /= 1) then
      write(*,*) "divide_class Wrong classes for C_1"
      stop
    end if
!
!  C_1
!
  ELSEIF (code_group==2.OR.code_group==3.OR.code_group==4) THEN
!
!  C_i, C_s, C_2
!
    IF (nclass /= 2) then
      write(*,*) "divide_class Wrong classes for C_i, C_s or C_2"
      stop
    end if
    which_irr(2)=2

  ELSEIF (code_group==5) THEN
!
!  C_3
!
! The function angle_rot(smat) provides the rotation angle of the matrix smat
!




end subroutine divide_class


subroutine type_symm(s,type_symm)
!--------------------------------------------------------------------------
! This function receives a 3x3 orthogonal matrix which is a symmetry
! operation of the point group of the crystal written in cartesian
! coordinates and gives as output a code according to the following:
!
!  1   Identity
!  2   Inversion
!  3   Proper rotation of an angle <> 180 degrees
!  4   Proper rotation of 180 degrees
!  5   Mirror symmetry
!  6   Improper rotation
!
  implicit none

  real(8), intent(in) :: s(3,3)
  real(8) :: det, det1,esp
  integer(4), intent(out) :: type_symm

  eps = EPS7

!
! Check for identity
!
  if ((abs(s(1,1)-1.d0) < eps) .and. &
      (abs(s(1,2) < eps) .and. &
      (abs(s(1,3) < eps) .and. &
      (abs(s(2,1) < eps) .and. &
      (abs(s(2,2)-1.d0) < eps) .and. &
      (abs(s(2,3) < eps) .and. &
      (abs(s(3,1) < eps) .and. &
      (abs(s(3,2) < eps) .and. &
      (abs(s(3,3)-1.d0) < eps)) then
   type_symm = 1
   return
  end if

!
! Check for inversion
!

  if ((abs(s(1,1)+1.d0) < eps) .and. &
      (abs(s(1,2) < eps) .and. &
      (abs(s(1,3) < eps) .and. &
      (abs(s(2,1) < eps) .and. &
      (abs(s(2,2)+1.d0) < eps) .and. &
      (abs(s(2,3) < eps) .and. &
      (abs(s(3,1) < eps) .and. &
      (abs(s(3,2) < eps) .and. &
      (abs(s(3,3)+1.d0) < eps)) then
    type_symm = 2
  end if

!
! compute the determinant
!
  det = s(1,1) * ( s(2,2) * s(3,3) - s(3,2) * s(2,3) )-   &
        s(1,2) * ( s(2,1) * s(3,3) - s(3,1) * s(2,3) )+   &
        s(1,3) * ( s(2,1) * s(3,2) - s(3,1) * s(2,2) )

  if (abs(det-1.d0) < eps) then
!
!  check if an eigenvalue is equal to -1.d0 (180 rotation)
!
   det1=(s(1,1)+1.d0)*((s(2,2)+1.d0)*(s(3,3)+1.d0)-s(3,2)*s(2,3))-   &
         s(1,2)*       (s(2,1)*      (s(3,3)+1.d0)-s(3,1)*s(2,3))+   &
         s(1,3)*       (s(2,1)*s(3,2)             -s(3,1)*(s(2,2)+1.d0))

   if (abs(det1) < eps) then
      type_symm = 4     ! 180 proper rotation
   else
      type_symm = 3     ! proper rotation <> 180
   endif
   return
  endif

!
! Determinant equal to -1: mirror symmetry or improper rotation
!
  if (abs(det+1.d0) < eps) then
!
!  check if an eigenvalue is equal to 1.d0 (mirror symmetry)
!
    det1=(s(1,1)-1.d0)*((s(2,2)-1.d0)*(s(3,3)-1.d0)-s(3,2)*s(2,3))-   &
          s(1,2)*       (s(2,1)*      (s(3,3)-1.d0)-s(3,1)*s(2,3))+   &
          s(1,3)*       (s(2,1)*s(3,2)             -s(3,1)*(s(2,2)-1.d0))

    if (abs(det1) < eps) then
      type_symm = 5   ! mirror symmetry
    else
      type_symm = 6   ! improper rotation
    endif
    return
  else
    write(*,*) "error symmetry not recognized."
    stop
  endif


end subroutine type_symm

end module symmetries
