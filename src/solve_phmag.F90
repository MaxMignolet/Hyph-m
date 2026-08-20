! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

module m_solve_phmag

USE m_commons
USE stdlib_sorting, ONLY: sort_index
USE m_lapack_wrappper
USE m_dynmat
USE m_phonons

IMPLICIT NONE

PRIVATE

PUBLIC :: solve_phmag

contains

subroutine solve_phmag(natom,nmagpert,amus,natom_species,atom_species,&
                      &fname_kpp,fname_kss,fname_kps,&
                      &fname_gpp,fname_gss,fname_gps)
!! FUNCTION
!! Compute coupled phonon-magnon frequencies
!!
!! INPUTS
!! natom: number of atoms
!! nmagpert: number of magnetic perturbations
!! amus: mass of the atoms, in amu
!! fname_***: name of the corresponding input files
!!
!! OUTPUT
!! Nothing
!!
!! SIDE EFFECT
!! Compute the coupled phonon-magnon frequencies

!Arguments ------------------------------------
integer,intent(in) :: natom,nmagpert,natom_species
real(dp),intent(in) :: amus(natom)
integer,intent(in) :: atom_species(natom)
CHARACTER(len=*) :: fname_kpp,fname_kss,fname_kps
CHARACTER(len=*) :: fname_gpp,fname_gss,fname_gps
!Local variables ------------------------------
! Scalar
integer :: io,i,ii,irow,imode,nat3,nat6,n
integer :: pmode,nmode,i_signed
CHARACTER(len=*),PARAMETER :: phmag_char_string(4) = [' P  ', '  M ', '   H', 'AP  ']
! Array
COMPLEX(dp),ALLOCATABLE :: k_pp(:,:) ! ph-ph
COMPLEX(dp),ALLOCATABLE :: g_pp(:,:)
COMPLEX(dp),ALLOCATABLE :: k_ss(:,:) ! sp-sp
COMPLEX(dp),ALLOCATABLE :: g_ss(:,:)
COMPLEX(dp),ALLOCATABLE :: k_ps(:,:) ! ph-sp
COMPLEX(dp),ALLOCATABLE :: g_ps(:,:)
COMPLEX(dp),ALLOCATABLE :: k_sp(:,:) ! sp-ph
COMPLEX(dp),ALLOCATABLE :: g_sp(:,:)
! Lapack Arrays
COMPLEX(dp),ALLOCATABLE :: vr(:,:)
COMPLEX(dp),ALLOCATABLE :: eigvec(:,:),eigvec_ucpld(:,:),eigvec_bare(:,:)
REAL(dp),ALLOCATABLE :: w(:),w_ucpld(:),w_bare(:)
! Other Arrays
REAL(dp),ALLOCATABLE :: phmag_character(:,:)
REAL(dp),ALLOCATABLE :: atomic_character(:,:)
REAL(dp),ALLOCATABLE :: phangmom_atomic(:,:,:)
real(dp) :: phmag_character_tot(3), phmag_character_mnus_tot(3), phmag_character_plus_tot(3)
REAL(dp),ALLOCATABLE :: phangmom(:,:)
! *************************************************************************

nat3=natom*3
nat6=natom*6
n=nat6+nmagpert

! read files
WRITE(*,*) "Read input files (once again)"


WRITE(*,*) "Reading k_pp..." ! Units: Ha/Bohr^2, but read as eV/Angstrom^2
ALLOCATE(k_pp(3*natom,3*natom))
open(newunit=io, file=fname_kpp)
do irow=1,3*natom
  read(io,*) k_pp(irow,1:3*natom)
  k_pp(irow,1:3*natom) = k_pp(irow,1:3*natom) / (27.2114/(0.52917)**2)
end do
call apply_asr_dynmat_q0(k_pp,natom)
call symmetrize_hermitian(k_pp,3*natom)
close(io)

WRITE(*,*) "Reading k_ps..." ! Units: Ha/Bohr, but read as eV/Angstrom
ALLOCATE(k_ps(3*natom,nmagpert))
open(newunit=io, file=fname_kps)
do irow=1,3*natom
  read(io,*) k_ps(irow,1:nmagpert)
  k_ps(irow,1:nmagpert) = k_ps(irow,1:nmagpert) / (27.2114/0.52917)
end do
close(io)
ALLOCATE(k_sp(nmagpert,nat3))
k_sp(:,:) = TRANSPOSE(k_ps)

WRITE(*,*) "Reading k_ss..." ! Units: Ha, but read as eV
ALLOCATE(k_ss(nmagpert,nmagpert))
open(newunit=io, file=fname_kss)
do irow=1,nmagpert
  read(io,*) k_ss(irow,1:nmagpert)
  k_ss(irow,1:nmagpert) = k_ss(irow,1:nmagpert) / 27.2114
end do
call symmetrize_hermitian(k_ss, nmagpert)
close(io)

WRITE(*,*) "Reading g_pp..." ! Units: 1/Bohr^2, but read as 1/Angstrom^2
ALLOCATE(g_pp(3*natom,3*natom))
open(newunit=io, file=fname_gpp)
do irow=1,3*natom
  read(io,*) g_pp(irow,1:3*natom)
  g_pp(irow,1:3*natom) = g_pp(irow,1:3*natom) / (1/(0.52917)**2)
end do
call symmetrize_antihermitian(g_pp,3*natom)
close(io)

WRITE(*,*) "Reading g_ps..." ! Units: 1/Bohr, but read as 1/Angstrom
ALLOCATE(g_ps(3*natom,nmagpert))
open(newunit=io, file=fname_gps)
do irow=1,3*natom
  read(io,*) g_ps(irow,1:nmagpert)
  g_ps(irow,1:nmagpert) = g_ps(irow,1:nmagpert) / (1/0.52917)
end do
close(io)
ALLOCATE(g_sp(nmagpert,nat3))
g_sp(:,:) = - TRANSPOSE(g_ps)

WRITE(*,*) "Reading g_ss..." ! Units: -
ALLOCATE(g_ss(nmagpert,nmagpert))
open(newunit=io, file=fname_gss)
do irow=1,nmagpert
  read(io,*) g_ss(irow,1:nmagpert)
end do
call symmetrize_antihermitian(g_ss, nmagpert)
close(io)

! multiplication by square root of the masses
call massmult(k_pp, natom, amus)
call massmult(g_pp, natom, amus)
call massmult_left(k_ps, natom, nmagpert, amus)
call massmult_left(g_ps, natom, nmagpert, amus)
call massmult_right(k_sp, nmagpert, natom, amus)
call massmult_right(g_sp, nmagpert, natom, amus)
! nothing to do for k_ss and g_ss

ALLOCATE(w(n))
ALLOCATE(w_ucpld(n))
ALLOCATE(w_bare(n))
ALLOCATE(eigvec(n,n))
ALLOCATE(eigvec_ucpld(n,n))
ALLOCATE(eigvec_bare(n,n))
ALLOCATE(vr(n,n))
ALLOCATE(phmag_character(3,n))
ALLOCATE(atomic_character(natom_species,n))
! Solving full system
call solve_system(natom,nmagpert,n,&
                 &k_pp,k_ss,k_ps,k_sp,&
                 &g_pp,g_ss,g_ps,g_sp,&
                 &w,eigvec,phmag_character,atomic_character,&
                 &atom_species)
! Solving uncoupled ph-mag
k_sp(:,:)=zero
k_ps(:,:)=zero
g_sp(:,:)=zero
g_ps(:,:)=zero
call solve_system(natom,nmagpert,n,&
                 &k_pp,k_ss,k_ps,k_sp,&
                 &g_pp,g_ss,g_ps,g_sp,&
                 &w_ucpld,eigvec_ucpld) ! _ucpld = uncoupled
! Solving bare uncoupled ph-mag
g_pp(:,:)=zero
call solve_system(natom,nmagpert,n,&
                 &k_pp,k_ss,k_ps,k_sp,&
                 &g_pp,g_ss,g_ps,g_sp,&
                 &w_bare,eigvec_bare)

WRITE(*,*) "Magnon-phonon frequencies | correction: bare (meV) | uncoupled (meV)"
do i=1,n/2
  imode=n/2+i
  WRITE(*,'(i3,2x,f14.8,f14.8,f14.8)')\
    i, w(imode)*Ha_2_meV,\
    (w(imode)-w_bare(imode))*Ha_2_meV,\
    (w(imode)-w_ucpld(imode))*Ha_2_meV
end do
WRITE(*,*) "Magnon-phonon frequencies: full |  bare | uncoupled (meV)"
do i=1,n/2
  imode=n/2+i
  WRITE(*,'(i3,2x,f14.8,f14.8,f14.8)')\
    i, w(imode)*Ha_2_meV,w_bare(imode)*Ha_2_meV,\
    w_ucpld(imode)*Ha_2_meV
end do
DEALLOCATE(w)
! Bare and ucpld quantities are not to be used again
DEALLOCATE(w_bare)
DEALLOCATE(w_ucpld)
DEALLOCATE(eigvec_bare)
DEALLOCATE(eigvec_ucpld)

WRITE(*,'(7x,a)') "Phononic    | Magnonic    | Coupled character (%):"
! do ii=n/2+1,n ! Only positive freq modes
do ii=1,n
  if (ii<=n/2) then
    i_signed = -(n/2-ii+1)
  else
    i_signed = ii-n/2
  endif
  WRITE(*,'(i3,2x,f14.8,f14.8,f14.8)') i_signed, phmag_character(:,ii)*100
end do

phmag_character_mnus_tot(:) = SUM(phmag_character(:,    1:n/2), DIM=2)
phmag_character_plus_tot(:) = SUM(phmag_character(:,n/2+1:n)  , DIM=2)
phmag_character_tot(:) = -phmag_character_mnus_tot(:) + phmag_character_plus_tot(:)
WRITE(*,'(a)') "Totals:"
WRITE(*,'(1x,a5,1x,3(f14.10),f15.10)') 'Mnus:', phmag_character_mnus_tot(:), SUM(phmag_character_mnus_tot)
WRITE(*,'(1x,a5,1x,3(f14.10),f15.10)') 'Plus:', phmag_character_plus_tot(:), SUM(phmag_character_plus_tot)
WRITE(*,'(1x,a5,1x,3(f14.10),f15.10)') 'Tot :', phmag_character_tot     (:), SUM(phmag_character_tot     )
WRITE(*,*) "Expected: ", nat6, nmagpert

! Computing/printing the phonon angular momentum
phangmom_block: block
  ! integer :: iat
  ! real(dp) :: pam_norm
  ! real(dp) :: pam_x,pam_y,pam_z,pam_xy
  ! real(dp) :: tilt_angle,roll_angle
  ALLOCATE(phangmom_atomic(3,natom,n)) ! cart dir, atom, mode
  call compute_phangmom_phmag(eigvec(1:nat3,:), eigvec(nat3+1:nat6,:), phangmom_atomic, natom, n)
  ALLOCATE(phangmom(3,n))
  phangmom(:,:) = sum(phangmom_atomic(:,:,:),dim=2) ! sum over atomic contributions

  WRITE(*,'(a)') "Phonon angular momentum (in units of hbar and in cart. coord):"
  WRITE(*,'(a)') "Printing pos.+neg. value"
  do imode=1,n/2
    pmode = n/2 + imode
    nmode = n/2 - imode + 1
    if (imode <= 3) then
      ii=4 ! Acoustic phonon
    else if (abs(phmag_character(1,pmode)) > half) then
      ii=1 ! phonon-like
    else if (abs(phmag_character(2,pmode)) > half) then
      ii=2 ! magnon-like
    else if (abs(phmag_character(3,pmode)) > half) then
      ii=3 ! hybrid-like
    endif
    WRITE(*,'(1x,i2,2x,3(es18.9,2x),a)')\
      imode,\
      phangmom(:,pmode) + phangmom(:,nmode),\
      trim(phmag_char_string(ii))
  end do
  WRITE(*,'(a,3(es18.9,2x))') "Net phonon angular momentum: ", sum(phangmom(:,:),dim=2)

  WRITE(*,'(a)') "Magnon angular momentum (in units of hbar and along z):"
  WRITE(*,'(a)') "Printing pos. mode values"
  do imode=1,n/2
    pmode = n/2 + imode
    nmode = n/2 - imode + 1
    if (imode <= 3) then
      ii=4 ! Acoustic phonon
    else if (abs(phmag_character(1,pmode)) > half) then
      ii=1 ! phonon-like
    else if (abs(phmag_character(2,pmode)) > half) then
      ii=2 ! magnon-like
    else if (abs(phmag_character(3,pmode)) > half) then
      ii=3 ! hybrid-like
    endif
    WRITE(*,'(1x,i2,2x,es18.9,2x,a)')\
      imode,\
      -phmag_character(2,pmode),\
      trim(phmag_char_string(ii))
  end do
  WRITE(*,'(a,es18.9,2x)') "Net magnon angular momentum: ", sum(phmag_character(2,1:n/2))

  ! ! Further post-processing of phangmom
  ! WRITE(*,'(a)') "Atomic resolved phangmom and tilt angle:"
  ! do ii=1,n/2
  !   imode = n/2 + ii
  !   if (norm2(phangmom(:,imode)) < tol2) cycle
  !   WRITE(*,*) "Mode: ", ii
  !   do iat=1,natom
  !     pam_norm = norm2(phangmom_atomic(:,iat,imode)*2)
  !     pam_x = phangmom_atomic(1,iat,imode)
  !     pam_y = phangmom_atomic(2,iat,imode)
  !     pam_z = phangmom_atomic(3,iat,imode)
  !     pam_xy = sqrt(pam_x**2 + pam_y**2)
  !     tilt_angle = atan2d(pam_xy,pam_z)
  !     roll_angle = atan2d(pam_y,pam_x)
  !     WRITE(*,'(1x,i2,2x,es18.9,2x,2(f14.8,2x))') iat, pam_norm, tilt_angle, roll_angle
  !   end do
  ! end do
end block phangmom_block

! eigvec_output_block: block
! integer :: iat
! do imode=1,n/2
!   WRITE(*,*) "Mode: ", imode
!   pmode = n/2 + imode
!   do iat=1,natom
!     WRITE(*,'(i2,2x,6(f12.8,2x))') iat, eigvec(3*(iat-1)+1:3*iat,pmode)
!   end do
! end do
! end block eigvec_output_block
DEALLOCATE(eigvec)

WRITE(*,'(a)') "Phonon atomic character (in % of total contribution, sum/=100%):"
WRITE(*,'(a)') "Printing pos. modes only"
WRITE(*,'(a)') " Mode   Cr              I/Br"
do imode=1,n/2

  pmode = n/2 + imode
  nmode = n/2 - imode + 1
  if (imode <= 3) then
    ii=4 ! Acoustic phonon
  else if (abs(phmag_character(1,pmode)) > half) then
    ii=1 ! phonon-like
  else if (abs(phmag_character(2,pmode)) > half) then
    ii=2 ! magnon-like
  else if (abs(phmag_character(3,pmode)) > half) then
    ii=3 ! hybrid-like
  endif
  WRITE(*,'(1x,i2,2x,2(f14.8,2x),a)')\
    imode,\
    atomic_character(:,n/2+imode)*100,\
    trim(phmag_char_string(ii))
end do

DEALLOCATE(phmag_character)
DEALLOCATE(atomic_character)
! Free all the k_xx
DEALLOCATE(k_pp)
DEALLOCATE(k_ss)
DEALLOCATE(k_ps)
DEALLOCATE(k_sp)
! Free all the g_xx
DEALLOCATE(g_pp)
DEALLOCATE(g_ss)
DEALLOCATE(g_ps)
DEALLOCATE(g_sp)

end subroutine solve_phmag

subroutine solve_system(natom,nmagpert,n,&
                       &k_pp,k_ss,k_ps,k_sp,&
                       &g_pp,g_ss,g_ps,g_sp,&
                       &w,eigvec,&
                       &phmag_character,atomic_character,&
                       &atom_species)
!! FUNCTION
!! Solve the system for coupled phonon-magnon
!!
!! INPUTS
!! natom: number of atoms
!! nmagpert: number of magnetic perturbations
!! n: natom*6+nmagpert ! it's just here for my convenience
!! k_**,g_**: IFC matrices and Brry curvature matrices
!!
!! OUTPUT
!! w: frequencies
!! eigvec: eigenvector of the system
!! phmag_character: for each mode gives the contribution of the purely atomic,
!! purely spin, and coupled atom-spin parts to the energy of the mode
!!
!! The frequencies and eigenvectors are sorted according to the frequencies
!! The eigenvectors are normalized to +\- i hbar
!! The ph/mag/coupled character is computed throught the M-norm
!!

!Arguments ------------------------------------
integer,intent(in) :: natom,nmagpert,n
COMPLEX(dp),intent(in) :: k_pp(natom*3,natom*3) ! ph-ph
COMPLEX(dp),intent(in) :: g_pp(natom*3,natom*3)
COMPLEX(dp),intent(in) :: k_ss(nmagpert,nmagpert) ! sp-sp
COMPLEX(dp),intent(in) :: g_ss(nmagpert,nmagpert)
COMPLEX(dp),intent(in) :: k_ps(natom*3,nmagpert) ! ph-sp
COMPLEX(dp),intent(in) :: g_ps(natom*3,nmagpert)
COMPLEX(dp),intent(in) :: k_sp(nmagpert,natom*3) ! sp-ph
COMPLEX(dp),intent(in) :: g_sp(nmagpert,natom*3)
REAL(dp),intent(out)    :: w(natom*6+nmagpert)
COMPLEX(dp),intent(out) :: eigvec(natom*6+nmagpert,natom*6+nmagpert)
REAL(dp),optional,intent(out) :: phmag_character(3,natom*6+nmagpert)
REAL(dp),optional,intent(out) :: atomic_character(:,:) ! (natom_species,natom*6+nmagpert)
integer,optional,intent(in)  :: atom_species(natom)

!Local variables ------------------------------
! Scalar
integer :: i,imode,nat3,nat6
! Array
integer :: ind(n)
! Dynamical matrices
COMPLEX(dp) :: D_pp(natom*3,natom*3)
COMPLEX(dp) :: D_ps(natom*3,nmagpert)
COMPLEX(dp) :: D_sp(nmagpert,natom*3)
COMPLEX(dp) :: D_ss(nmagpert,nmagpert)
! Lapack Arrays
COMPLEX(dp) :: matA(n,n)
COMPLEX(dp) :: matB(n,n), LHS(n,n) ! LHS will be a copy of matB
COMPLEX(dp) :: alpha(n)
COMPLEX(dp) :: beta(n)
COMPLEX(dp) :: vr(n,n)
! *************************************************************************

nat3=natom*3
nat6=natom*6

! Forming the dynamical matrices
D_pp(:,:) = k_pp - matmul(g_pp,g_pp) /4
D_ps(:,:) = k_ps - matmul(g_pp,g_ps) /4
D_sp(:,:) = k_sp - matmul(g_sp,g_pp) /4
D_ss(:,:) = k_ss - matmul(g_sp,g_ps) /4

! Solving the equations now...
! -iw (  0    -1    G_us) (u) = ( D_pp,    -G_pp,   D_ps)  (u)
! -iw (  1     0    0   ) (p) = ( G_pp,     1,      G_ps)  (p)
! -iw (  G_su  0    G_ss) (s) = ( D_sp,    -G_sp,   D_ss)  (s)
! All the G matrices are implicitely divided by two except G_ss
! left hand side correspond the lapack B matrix
! right hand side correspond the lapack A matrix
WRITE(*,*) "Construction of the big matrices..."
matA(:,:) = zero
matB(:,:) = zero

! Left hand side first
! First row
do i=1,nat3
  matB(i,nat3+i) = -cone
end do
matB(1:nat3,nat6+1:n) = g_ps /2
! Second row
do i=1,nat3
  matB(nat3+i,i) = cone
end do
! Third row
matB(nat6+1:n,1:nat3) = g_sp /2
matB(nat6+1:n,nat6+1:n) = g_ss

! Right hand side
! First row
matA(1:nat3,     1:nat3) =  D_pp
matA(1:nat3,nat3+1:nat6) = -g_pp /2
matA(1:nat3,nat6+1:n   ) =  D_ps
! Second row
matA(nat3+1:nat6,1:nat3) = g_pp /2
do i=1,nat3
  matA(nat3+i,nat3+i) = cone
end do
matA(nat3+1:nat6,nat6+1:n) = g_ps /2
! Third row
matA(nat6+1:n,     1:nat3) = D_sp
matA(nat6+1:n,nat3+1:nat6) = -g_sp /2
matA(nat6+1:n,nat6+1:n) = D_ss

WRITE(*,'(a)') "Solving system..."
LHS(:,:) = matB(:,:) ! make a copy before the call (will be used for the
                     ! normalization of the eigvec)
call zggev3_vr_wrp(matA, matB, n, alpha, beta, vr)
! For some reason I can't use eigvec directly, I have to go through another array
WRITE(*,'(a)') "Solved."

w(:) = -AIMAG(alpha/beta)
! Sorting modes according to their frequencies
! The zero frequency modes are not necessarily correctly sorted..
call sort_index(w, ind)

! Copying and sorting the eigenvectors
eigvec(:,:) = vr(:,:)
do imode=1,n
  vr(:,imode) = eigvec(:,ind(imode))
end do
eigvec(:,:) = vr(:,:)


! computation of the contribution of the ph/magnon/coupled parts to the mode
! energy for each mode
!                 ( 0    -1    G_us) (u)
! (u^H  p^H  s^H) ( 1     0    0   ) (p) = +/-i hbar
!                 ( G_su  0    G_ss) (s)
!                (0  -1) (u)
! ph: (u^H  p^H) (1   0) (p)
!
! mag: s^H G_ss s
!                    (0    G_us) (u)
! coupled: (u^H s^H) (G_su 0   ) (s)
!
call renorm_phmag(eigvec,LHS,n)

! if we don't have to compute the phmag_character, we can stop here and return
if (.not. present(phmag_character)) return
if (.not. present(atomic_character)) return
if (.not. present(atom_species)) return

! Compute phonon-magnon character
call compute_phmag_charac_per_mode(eigvec,natom,nmagpert,n,atom_species,LHS,&
                                   phmag_character,atomic_character)

end subroutine solve_system

end module m_solve_phmag
