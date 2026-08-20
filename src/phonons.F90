! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

MODULE m_phonons

use m_commons
use stdlib_sorting, only: sort_index

implicit none


private

public :: sort_modes
public :: renorm_phmbc
public :: renorm_phmag_Ren24
public :: renorm_phmag
public :: compute_phangmom_phmbc
public :: compute_phangmom_phmag
public :: compute_atomic_charac_per_mode
public :: compute_atomic_charac_per_phmbc_mode
public :: compute_phmag_charac_per_mode_SRen24
public :: compute_phmag_charac_per_mode
public :: phdispl_from_pheigvec
public :: gen_ph_trajectory
public :: dump_ph_trajectory

CONTAINS


subroutine sort_modes(freq,eigvec,n,nmode)
!! FUNCTION
!! Sort modes according to their frequency
!!
!! INPUTS
!! freq(nmode): frequenciess fo the mode
!! eigvec(n,nmode): eigenvector of the corresponding mode
!! n: size of the eigenvectors
!! nmode: number of modes
!!
!! OUTPUT
!! freq and eigvec sorted

!Arguments ------------------------------------
integer,intent(in) :: n,nmode
real(dp),intent(inout) :: freq  (nmode)
complex(dp),intent(inout) :: eigvec(n,nmode)
!Local variables ------------------------------
integer :: imode
integer :: ind(n)
complex(dp) :: tmp(n,nmode)
! *************************************************************************

! Sorting modes according to their frequencies
! as a word of caution the zero frequency modes are not necessarily correctly
! sorted. If you need to reuse them, be careful..
call sort_index(freq, ind)
freq(:) = freq(ind(:))

! Sorting the eigenvectors
do imode=1,n
  tmp(:,imode) = eigvec(:,ind(imode))
end do
eigvec(:,:) = tmp(:,:)

end subroutine sort_modes
!!***

subroutine renorm_phmbc(pheigvec,pheigvec_vel,natom,nmode)
!! FUNCTION
!! Normalize phonons modes according to
!! (pheigvec^* . pheigvec_vel + pheigvec_vel^* . pheigvec)/2 = 1
!!
!! INPUTS
!! pheigvec: phonon displacement eigenvector in cart coord
!! pheigvec_vel: phonon velocity eigenvector in cart coord
!! natom: number of atoms
!! nmode: number of modes
!!
!! OUTPUT
!! pheigvec and pheivec_vel renormalized
!!
!! Notes
!! Used in the ph+mbc system

!Arguments ------------------------------------
integer,intent(in) :: natom,nmode
complex(dp),intent(inout) :: pheigvec    (3*natom,nmode)
complex(dp),intent(inout) :: pheigvec_vel(3*natom,nmode)
!Local variables ------------------------------
integer :: imode
real(dp) :: norm2, norm_displ
! *************************************************************************

! normalize the 3 acoustic modes for which the norm would be zero
do imode=1,3
  ! normalize displacement part to one, and set the velocity to zero
  norm_displ = REAL(DOT_PRODUCT(pheigvec(:,imode), pheigvec(:,imode)))
  pheigvec    (:,imode) = pheigvec(:,imode) / sqrt(norm_displ)
  pheigvec_vel(:,imode) = 0
end do

do imode=4,nmode
  norm2 = REAL(DOT_PRODUCT(pheigvec(:,imode), pheigvec_vel(:,imode))\
        + DOT_PRODUCT(pheigvec_vel(:,imode), pheigvec(:,imode)),dp) / 2
  pheigvec    (:,imode) = pheigvec    (:,imode) / sqrt(norm2)
  pheigvec_vel(:,imode) = pheigvec_vel(:,imode) / sqrt(norm2)

  ! normalize displacement part
  norm_displ = REAL(DOT_PRODUCT(pheigvec(:,imode), pheigvec(:,imode)))
  pheigvec    (:,imode) = pheigvec    (:,imode) / sqrt(norm_displ)
  pheigvec_vel(:,imode) = pheigvec_vel(:,imode) * sqrt(norm_displ)
end do

end subroutine renorm_phmbc
!!***

subroutine renorm_phmag_Ren24(eigvec,natom,n)
!! FUNCTION
!! Normalize phonon-magnon modes according to
!! |phdispl|**2 + |speigvec|**2 = 1
!!
!! INPUTS
!! eigvec: eigenvector in cart coord
!! natom: number of atoms
!! nmagpert: number of magnetic perturbations
!! n: natom*6 + nmagpert
!!
!! OUTPUT
!! eigvec renormalized
!!
!! Notes
!! Used in the 1st method for ph-mag

!Arguments ------------------------------------
integer,intent(in) :: natom,n
complex(dp),intent(inout) :: eigvec(n,n)
!Local variables ------------------------------
integer :: imode
complex(dp) :: norm_cplx
real(dp) :: norm_mode
! *************************************************************************

do imode=1,n
  norm_cplx = DOT_PRODUCT(eigvec(1:natom*3,imode),   eigvec(1:natom*3,imode))&
            + DOT_PRODUCT(eigvec(natom*6+1:n,imode), eigvec(natom*6+1:n,imode))
  norm_mode = real(norm_cplx)

  ! the six modes with the smallest freq are the acoustic ones
  ! i.e. the six mode in the middle of the spectrum
  if (n/2-3 < imode .and. imode <= n/2+3) then
    cycle
  endif

  eigvec(:,imode) = eigvec(:,imode) / sqrt(abs(norm_mode))
end do

end subroutine renorm_phmag_Ren24
!!***

subroutine renorm_phmag(eigvec,M,n)
!! FUNCTION
!! Normalize phonon-magnon modes according to
!! eigvec^* M eigvec = +/- i hbar
!!
!! INPUTS
!! eigvec: eigenvector in cart coord
!! M matrix: left hand side
!! natom: number of atoms
!! nmagpert: number of magnetic perturbations
!! n: natom*6 + nmagpert
!!
!! OUTPUT
!! eigvec renormalized
!!
!! Notes
!! Used in the full ph-mag system

!Arguments ------------------------------------
integer,intent(in) :: n
complex(dp),intent(inout) :: eigvec(n,n)
complex(dp),intent(in) :: M(n,n)
!Local variables ------------------------------
integer :: imode
complex(dp) :: norm_cplx
real(dp) :: norm_mode
! *************************************************************************

do imode=1,n
  norm_cplx = DOT_PRODUCT(eigvec(:,imode), MATMUL(M, eigvec(:,imode))) / j_dpc
  norm_mode = real(norm_cplx)

  ! the six modes with the smallest freq are the acoustic ones
  ! i.e. the six mode in the middle of the spectrum
  if (n/2-3 < imode .and. imode <= n/2+3) then
    cycle
  endif

  eigvec(:,imode) = eigvec(:,imode) / sqrt(abs(norm_mode))

  ! norm_cplx = DOT_PRODUCT(eigvec(1:3*nat,imode),eigvec(1:3*nat,imode))
  ! norm_displ = real(norm_cplx)
  ! eigvec(1:3*nat,imode) = eigvec(1:3*nat,imode) / sqrt(abs(norm_displ))
  ! eigvec(3*nat+1:6*nat,imode) = eigvec(3*nat+1:6*nat,imode) * sqrt(abs(norm_displ))
end do

end subroutine renorm_phmag
!!***

subroutine compute_phangmom_phmbc(pheigvec,pheigvec_vel,phangmom,natom,nmode)
!! FUNCTION
!! Compute the angular momentum of phonon modes. The modes should be normalized
!! via renorm_phonon
!! L_{nu,kappa} = Im{pheigvec^* x pheigvec_vel}
!!
!! INPUTS
!! pheigvec: phonon displacement eigenvector in cart coord
!! pheigvec_vel: phonon velocity eigenvector in cart coord
!! natom: number of atoms
!! nmode: number of modes
!!
!! OUTPUT
!! phangmom: contains the angular momentum

!Arguments ------------------------------------
integer,intent(in) :: natom,nmode
complex(dp),intent(in) :: pheigvec    (3,natom,nmode)
complex(dp),intent(in) :: pheigvec_vel(3,natom,nmode)
real(dp),intent(out) :: phangmom(3,natom,nmode)
!Local variables ------------------------------
integer :: imode,iat
complex(dp) :: eigvec_atom(3)
complex(dp) :: eigvec_vel_atom(3)
! *************************************************************************

phangmom(:,:,:) = zero
do imode=1,nmode
  do iat=1,natom
    eigvec_atom(:)     = pheigvec    (:,iat,imode)
    eigvec_vel_atom(:) = pheigvec_vel(:,iat,imode)
    ! l_x = u_y p_z - u_z p_y
    phangmom(1,iat,imode) = aimag(conjg(eigvec_atom(2)) * eigvec_vel_atom(3)\
                                - conjg(eigvec_atom(3)) * eigvec_vel_atom(2))
    ! l_y = u_z p_x - u_x p_z
    phangmom(2,iat,imode) = aimag(conjg(eigvec_atom(3)) * eigvec_vel_atom(1)\
                                - conjg(eigvec_atom(1)) * eigvec_vel_atom(3))
    ! l_z = u_x p_y - u_y p_x
    phangmom(3,iat,imode) = aimag(conjg(eigvec_atom(1)) * eigvec_vel_atom(2)\
                                - conjg(eigvec_atom(2)) * eigvec_vel_atom(1))
  end do
end do

end subroutine compute_phangmom_phmbc
!!***

subroutine compute_phangmom_phmag(pheigvec,pheigvec_vel,phangmom,natom,nmode)
!! FUNCTION
!! Compute the angular momentum of phonon modes. The modes should be normalized
!! via renorm_phonon
!! L_{nu,kappa} = Re{ pheigvec^* x pheigvec_vel }
!!
!! INPUTS
!! pheigvec: phonon displacement eigenvector in cart coord
!! pheigvec_vel: phonon velocity eigenvector in cart coord
!! natom: number of atoms
!! nmode: number of modes
!!
!! OUTPUT
!! phangmom: contains the angular momentum

!Arguments ------------------------------------
integer,intent(in) :: natom,nmode
complex(dp),intent(in) :: pheigvec    (3,natom,nmode)
complex(dp),intent(in) :: pheigvec_vel(3,natom,nmode)
real(dp),intent(out) :: phangmom(3,natom,nmode)
!Local variables ------------------------------
integer :: imode,iat
complex(dp) :: eigvec_atom(3)
complex(dp) :: eigvec_vel_atom(3)
! *************************************************************************

phangmom(:,:,:) = zero
do imode=1,nmode
  do iat=1,natom
    eigvec_atom(:)     = pheigvec    (:,iat,imode)
    eigvec_vel_atom(:) = pheigvec_vel(:,iat,imode)
    ! l_x = u_y p_z - u_z p_y
    phangmom(1,iat,imode) = real(conjg(eigvec_atom(2)) * eigvec_vel_atom(3)\
                               - conjg(eigvec_atom(3)) * eigvec_vel_atom(2))
    ! l_y = u_z p_x - u_x p_z
    phangmom(2,iat,imode) = real(conjg(eigvec_atom(3)) * eigvec_vel_atom(1)\
                               - conjg(eigvec_atom(1)) * eigvec_vel_atom(3))
    ! l_z = u_x p_y - u_y p_x
    phangmom(3,iat,imode) = real(conjg(eigvec_atom(1)) * eigvec_vel_atom(2)\
                               - conjg(eigvec_atom(2)) * eigvec_vel_atom(1))
  end do
end do

end subroutine compute_phangmom_phmag
!!***

subroutine compute_atomic_charac_per_mode(pheigvec,natom,atom_species,natom_species,atomic_character)
!! FUNCTION
!! Compute the atomic character of each phonon mode, based on the phonon
!! displacmenet eigenvector in cartesian coordinates.
!!
!! INPUTS
!! pheigvec: phonon displacement vector in cart coord
!! natom: number of atoms
!! atom_species(natom_species): type of each atom
!! natom_species: number of different types
!!
!! OUTPUT
!! atomic_character(natom_species,3*natom): contribution of each atomic species
!! to each mode. The sum of the contributionns for one mode sums up to one

!Arguments ------------------------------------
integer,intent(in) :: natom,natom_species
integer,intent(in) :: atom_species(natom)
complex(dp),intent(in) :: pheigvec(3,natom,3*natom)
real(dp),intent(out) :: atomic_character(natom_species,3*natom)
!Local variables ------------------------------
integer :: imode,iat
! *************************************************************************

atomic_character = zero
do imode=1,3*natom
  do iat=1,natom
    atomic_character(atom_species(iat),imode) = \
      atomic_character(atom_species(iat),imode) + \
      REAL(DOT_PRODUCT(pheigvec(:,iat,imode), pheigvec(:,iat,imode)),dp)
  end do
end do

end subroutine compute_atomic_charac_per_mode
!!***

subroutine compute_atomic_charac_per_phmbc_mode(pheigvec,pheigvec_vel,natom,n,&
                                   &atom_species,natom_species,atomic_character)
!! FUNCTION
!! Compute the atomic character of each phonon+mbc mode, based on the phonon
!! eigenvector in cartesian coordinates.
!!
!! INPUTS
!! pheigvec: phonon vector in cart coord
!! natom: number of atoms
!! n: number of modes
!! atom_species(natom_species): type of each atom
!! natom_species: number of different types
!!
!! OUTPUT
!! atomic_character(natom_species,3*natom): contribution of each atomic species
!! to each mode. The sum of the contributionns for one mode sums up to one

!Arguments ------------------------------------
integer,intent(in) :: natom,n,natom_species
integer,intent(in) :: atom_species(natom)
complex(dp),intent(in) :: pheigvec    (3,natom,n)
complex(dp),intent(in) :: pheigvec_vel(3,natom,n)
real(dp),intent(out) :: atomic_character(natom_species,n)
!Local variables ------------------------------
integer :: imode,iat
complex(dp) :: displ(3), vel(3)
! *************************************************************************

atomic_character = zero
do imode=1,n
  do iat=1,natom
    displ(:) = pheigvec    (:,iat,imode)
    vel(:)   = pheigvec_vel(:,iat,imode)
    atomic_character(atom_species(iat),imode) = \
      atomic_character(atom_species(iat),imode) + \
      REAL(DOT_PRODUCT(displ, vel) + DOT_PRODUCT(vel, displ),dp)/2
  end do
end do

end subroutine compute_atomic_charac_per_phmbc_mode
!!***

subroutine compute_phmag_charac_per_mode_SRen24(eigvec,natom,nmagpert,n,phmag_character)
!! FUNCTION
!! Compute the phononic/magnonic character of each phonon-magnon mode, based on
!! the full phonon-magnon eigenvector in cartesian coordinates.
!! This does not 100% match what is in the paper, here I still take into account
!! the effect of the mbc on the normalization (whereas it is not done in the
!! paper, cfr Eq. 13 of their paper)
!!
!! INPUTS
!! eigvec: phonon-magnon vector in cart coord
!! natom: number of atoms
!! nmagpert: number of magnetic perturbations
!! n: natom*6 + natom
!!
!! OUTPUT
!! phmag_character(2,n): contribution of the phonon/magnon part to each
!! mode. The sum of the contributionns for one mode sums up to one

!Arguments ------------------------------------
integer,intent(in) :: natom,nmagpert,n
complex(dp),intent(in) :: eigvec(n,n)
real(dp),intent(out) :: phmag_character(2,n)
!Local variables ------------------------------
integer :: imode,nat3,nat6,norm_sign
complex(dp) :: pheigvec    (natom*3)
complex(dp) :: pheigvec_vel(natom*3)
complex(dp) :: speigvec    (nmagpert)
! *************************************************************************

nat3=3*natom
nat6=6*natom

phmag_character = zero
do imode=1,n
  if (imode<=n/2) then
    norm_sign = -1
  else
    norm_sign = 1
  end if
  pheigvec(:) = eigvec(1:nat3,imode)
  pheigvec_vel(:) = eigvec(nat3+1:nat6,imode)
  speigvec(:) = eigvec(nat6+1:n,imode)

  phmag_character(1,imode) = norm_sign * \
                             REAL(DOT_PRODUCT(pheigvec, pheigvec_vel)\
                                + DOT_PRODUCT(pheigvec_vel, pheigvec),dp) / 2
  phmag_character(2,imode) = REAL(DOT_PRODUCT(speigvec, speigvec),dp)
end do

end subroutine compute_phmag_charac_per_mode_SRen24
!!***

subroutine compute_phmag_charac_per_mode(eigvec,natom,nmagpert,n,atom_species,M_matrix,phmag_character,atomic_character)
!! FUNCTION
!! Compute the phononic/magnonic character of each phonon-magnon mode, based on
!! the full phonon-magnon eigenvector in cartesian coordinates.
!! Also compute the atomic character for each mode
!!
!! INPUTS
!! eigvec: phonon-magnon vector in cart coord
!! natom: number of atoms
!! nmagpert: number of magnetic perturbations
!! n: natom*6 + nmagpert
!!
!! OUTPUT
!! phmag_character(3,n): contribution of the phonon/magnon/cpld part to each
!! mode. The sum of the contributionns for one mode sums up to one

!Arguments ------------------------------------
integer,intent(in) :: natom,nmagpert,n
complex(dp),intent(in) :: eigvec(n,n)
integer,intent(in)  :: atom_species(natom)
COMPLEX(dp),intent(in) :: M_matrix(n,n)
real(dp),intent(out) :: phmag_character(3,n)
REAL(dp),optional,intent(out) :: atomic_character(:,:) ! (natom_species,natom*6+nmagpert)
!Local variables ------------------------------
integer :: imode,iat,nat3,nat6
complex(dp),allocatable :: phdispl(:), phvel(:), phvec(:), magdispl(:)
complex(dp) :: displ(3), vel(3)
real(dp) :: ph_char, mag_char, cpld_char
complex(dp) :: norm_cplx
! *************************************************************************

nat3=3*natom
nat6=6*natom
ALLOCATE(phdispl(nat3))
ALLOCATE(phvel(nat3))
ALLOCATE(phvec(nat6))
ALLOCATE(magdispl(nmagpert))

phmag_character = zero
atomic_character = zero
do imode=1,n
  ! the six modes with the smallest freq are the acoustic ones
  ! i.e. the six mode in the middle of the spectrum
  ! differentiated for neg/pos modes
  if (n/2-3 < imode .and. imode <= n/2) then
    phmag_character(1,imode) =-one  ! ph_char
    phmag_character(2,imode) = zero ! mag_char
    phmag_character(3,imode) = zero ! cpld_char
    cycle
  endif
  if (n/2   < imode .and. imode <= n/2+3) then
    phmag_character(1,imode) = one  ! ph_char
    phmag_character(2,imode) = zero ! mag_char
    phmag_character(3,imode) = zero ! cpld_char
    cycle
  endif
  phdispl (:) = eigvec(     1:nat3,imode)
  phvel   (:) = eigvec(nat3+1:nat6,imode)
  phvec   (:) = eigvec(1:nat6,imode)
  magdispl(:) = eigvec(nat6+1:n   ,imode)

  ! phonon contribution to energy
  norm_cplx = DOT_PRODUCT(phvec, MATMUL(M_matrix(1:nat6,1:nat6), phvec)) / j_dpc
  ph_char = real(norm_cplx)

  ! magnon contribution to energy
  norm_cplx = DOT_PRODUCT(magdispl, MATMUL(M_matrix(nat6+1:n,nat6+1:n), magdispl)) / j_dpc
  mag_char = real(norm_cplx)

  ! coupled contribution to energy
  norm_cplx = (DOT_PRODUCT(phdispl , MATMUL(M_matrix(1:nat3,nat6+1:n), magdispl))\
              +DOT_PRODUCT(magdispl, MATMUL(M_matrix(nat6+1:n,1:nat3), phdispl ))) / j_dpc
  cpld_char = real(norm_cplx)

  phmag_character(1,imode) = ph_char
  phmag_character(2,imode) = mag_char
  phmag_character(3,imode) = cpld_char

  do iat=1,natom
    displ(:) = phdispl(3*(iat-1)+1:3*iat)
    vel(:)   = phvel  (3*(iat-1)+1:3*iat)
    atomic_character(atom_species(iat),imode) = &
      atomic_character(atom_species(iat),imode) &
      + REAL((-DOT_PRODUCT(displ, vel) + DOT_PRODUCT(vel, displ))/j_dpc,dp)
  end do
end do

end subroutine compute_phmag_charac_per_mode
!!***

subroutine phdispl_from_pheigvec(natom, amus, pheigvec, phdispl)
!! The approach taken here seems outdated, the normalization irks me a bit
!! FUNCTION
!!  Phonon displacements in cart coords from eigenvectors for one mode
!!
!! INPUTS
!!  natom: number of atoms in unit cell
!!  amus(natom)=mass of the atoms (atomic mass unit)
!!  pheigvec(3,natom)= phonon displacement eigenvectors in cartesian coordinates.
!!
!! OUTPUT
!!  phdispl(3,natom)=displacements of atoms in cartesian coordinates.

!Arguments -------------------------------
!scalars
integer,intent(in) :: natom
!arrays
real(dp),intent(in) :: amus(natom)
complex(dp),intent(in) :: pheigvec(3,natom)
complex(dp),intent(out) :: phdispl(3,natom)
!Local variables -------------------------
!scalars
integer :: iat
real(dp) :: norm
! *********************************************************************

norm=0
do iat=1,natom
  phdispl(:,iat)=pheigvec(:,iat) / sqrt(amus(iat)*amu_emass)
  norm = norm + REAL(DOT_PRODUCT(phdispl(:,iat),\
                                 phdispl(:,iat)),dp)
end do
! renormalization
norm = sqrt(norm)
phdispl(:,:) = phdispl(:,:) / norm

end subroutine phdispl_from_pheigvec
!!***

subroutine gen_ph_trajectory(pheigvec,phfreq,xcart,natom,amus,&
                             timestep,nperiod,traj)
!! FUNCTION
!! Generate atomic trajectories based on the phonon displacmenet eigenvector in
!! cartesian coordinates. (for one specific mode)
!! This is intended to be used for visualization, and especially to be fed into
!! dump_ph_trajectory. The routine will compute the positions of the atoms
!! at `timestep` interval for `nperiod` periods (recommended values: 2-5)
!!
!! INPUTS
!! pheigvec: phonon displacement vector in cart coord for the wanted mode
!! phfreq: phonon frequency in Ha
!! natom: number of atoms
!! timestep: timestep in 1e-12 seconds
!! nperiod: number of periods to be generated
!!
!! OUTPUT
!! trajectory(natom,3*natom): trajectory of the atoms, it is a list of the
!! atomic positions in cart coord.

!Arguments ------------------------------------
integer,intent(in) :: natom,nperiod
real(dp),intent(in) :: timestep,xcart(3,natom),amus(natom)
real(dp),intent(in) :: phfreq
complex(dp),intent(in) :: pheigvec(3,natom)

real(dp),ALLOCATABLE,intent(out) :: traj(:,:,:)
!Local variables ------------------------------
integer :: nframe,iframe,iat
real(dp) :: tstep,period,time,displ_cart(3,natom),pos_cart(3,natom)
complex(dp) :: phdispl(3,natom)
! *************************************************************************


period = two_pi/(phfreq*Ha_THz)
nframe = int(nperiod * period/timestep) ! we want it to loop perfectly
tstep = nperiod*period/nframe ! new adapted timestep
time = 0

ALLOCATE(traj(3,natom,nframe))
call phdispl_from_pheigvec(natom, amus, pheigvec, phdispl)

do iframe=1,nframe
  time = (iframe-1)*timestep
  do iat=1,natom
    displ_cart(:,iat) = real(phdispl(:,iat) * exp(-j_dpc*phfreq*Ha_THz*time),dp)
  end do
  pos_cart(:,:) = xcart + displ_cart*0.8
  do iat=1,natom
    traj(:,iat,iframe) = pos_cart(:,iat)
  end do
end do

end subroutine gen_ph_trajectory
!!***

subroutine dump_ph_trajectory(prefix,latt_vec,natom,atom_species,natom_species,&
                              atom_species_label,traj)
!! FUNCTION
!! Dump atomic trajectories based on the phonon displacmenet eigenvector in
!! cartesian coordinates. (for one specific mode)
!! Trajectories should have been generated by gen_ph_trajectory
!!
!! INPUTS
!! latt_vec: lattice vector in Bohr
!! natom: number of atoms
!! atom_species(natom_species): type of each atom
!! natom_species: number of different types
!! traj(3,natom,nframe): atomic trajectories (see gen_ph_trajectory)
!!
!! OUTPUT
!! Writing to file

!Arguments ------------------------------------
character(len=*),intent(in) :: prefix
integer,intent(in) :: natom,natom_species
integer,intent(in) :: atom_species(natom)
real(dp),intent(in) :: latt_vec(3,3)
character(len=2),intent(in) :: atom_species_label(natom_species)
real(dp),intent(in) :: traj(:,:,:)
!Local variables ------------------------------
integer :: nframe,iframe,iat,unit
character(len=100) :: msg,fname
! *************************************************************************

nframe = SIZE(traj,dim=3)
do iframe=1,nframe
  WRITE(fname,'(a,i4.4)') trim(prefix), iframe
  open(newunit=unit,file=fname,action='write')

  ! Comment line
  WRITE(unit,'(a)') "CrI3"
  ! Lattice vectors
  WRITE(unit, '(es22.14,a)') Bohr_Angstrom
  WRITE(unit, '(3(es22.14,2x))') latt_vec(:,1)
  WRITE(unit, '(3(es22.14,2x))') latt_vec(:,2)
  WRITE(unit, '(3(es22.14,2x))') latt_vec(:,3)

  ! Write labels
  msg = ''
  do iat=1,natom_species
    WRITE(msg,'(a,2x,a2)') trim(msg), atom_species_label(iat)
  end do
  WRITE(unit, '(a)') trim(adjustl(msg))

  ! Write number of each atomic species
  msg = ''
  do iat=1,natom_species
    WRITE(msg,'(a,2x,i4)') trim(msg), count(atom_species==iat)
  end do
  WRITE(unit, '(a)') trim(adjustl(msg))

  WRITE(unit, '(a)') 'Cartesian'

  ! Write atomic positions
  do iat=1,natom
    WRITE(unit,'(3(es22.14,2x))') traj(:,iat,iframe)
  end do

  close(unit)
end do

end subroutine dump_ph_trajectory
!!***

END MODULE
