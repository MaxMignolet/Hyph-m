! Hyph-m - main routine

! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

PROGRAM MAIN

USE stdlib_sorting, ONLY: sort_index
USE m_commons
USE m_dynmat
USE m_lapack_wrappper
USE m_phonons
USE m_readInput
USE m_solve_saparov
USE m_solve_phmag

IMPLICIT NONE

! Variables
!! Parameters
!!! Phonon files
CHARACTER(len=*),PARAMETER :: fname_kpp="indata/k_pp.txt"
CHARACTER(len=*),PARAMETER :: fname_gpp="indata/g_pp.txt"
!!! Spin files
CHARACTER(len=*),PARAMETER :: fname_kss="indata/k_ss.txt"
CHARACTER(len=*),PARAMETER :: fname_gss="indata/g_ss.txt"
!!! Spin-Phonon coupling files
CHARACTER(len=*),PARAMETER :: fname_kps="indata/k_ps.txt"
CHARACTER(len=*),PARAMETER :: fname_gps="indata/g_ps.txt"
!!! Other parameters
INTEGER,PARAMETER :: prt_debug=0, max_prt_lines=40
CHARACTER(LEN=*),PARAMETER :: cart(3)=(/'x','y','z'/)
!! General IO
CHARACTER(LEN=len100) :: msg
INTEGER :: io

!! loops variables
INTEGER :: i,ii,jj,iat,imode,irow
INTEGER :: iat1,iat2
INTEGER :: idir1,idir2
INTEGER :: ipc1,ipc2
INTEGER,ALLOCATABLE :: ind(:)
! REAL(dp) :: norm

!! Atom species and masses
character(len=len100) :: crystName
INTEGER :: natom, nat3, nat6, natom_species
INTEGER,ALLOCATABLE :: atom_species(:)! natom, stores to species of each atom
real(dp),ALLOCATABLE :: atomic_masses(:) ! stores the mass of each atom species
character(len=2),ALLOCATABLE :: atom_species_label(:) ! label for each atomic species
INTEGER :: nmagatom ! number of magnetic atom
INTEGER,ALLOCATABLE :: magatom(:) ! their index
INTEGER :: nmagdir ! and the number of direction considered (here only x and y, so iit will be 2)
INTEGER :: nmagpert ! number of atomic magnetic perturbations (nmagatom*nmagdir)
real(dp),ALLOCATABLE :: xred(:,:)! atomic positions in red coord
real(dp),ALLOCATABLE :: xcart(:,:)
real(dp) :: latt_vec(3,3) ! lattice vectors in Angstrom
real(dp),ALLOCATABLE :: amus(:) ! natom, stores the mass of each atom

!! Eigenvalue problem variables
COMPLEX(dp),ALLOCATABLE :: k_pp(:,:) ! IFC matrix
COMPLEX(dp),ALLOCATABLE :: g_pp(:,:) ! MBC matrix
COMPLEX(dp),ALLOCATABLE :: k_ss(:,:) ! k_ss matrix
COMPLEX(dp),ALLOCATABLE :: g_ss(:,:) ! SBC matrix
COMPLEX(dp),ALLOCATABLE :: k_ps(:,:) !
COMPLEX(dp),ALLOCATABLE :: g_ps(:,:) !
COMPLEX(dp),ALLOCATABLE :: k_sp(:,:) !
COMPLEX(dp),ALLOCATABLE :: g_sp(:,:) !
COMPLEX(dp),ALLOCATABLE :: tmp_matrix(:,:)
REAL(dp),ALLOCATABLE :: phfreq(:) ! phonon frequencies
COMPLEX(dp),ALLOCATABLE :: eigvec(:,:) ! general eignvec
COMPLEX(dp),ALLOCATABLE :: pheigvec(:,:,:) ! phonon eigenvec(idir,iat,imode)
COMPLEX(dp),ALLOCATABLE :: pheigvec_vel(:,:,:) ! velocity part
REAL(dp),ALLOCATABLE :: phangmom(:,:) ! phonon angular momentum per mode
REAL(dp),ALLOCATABLE :: phangmom_atomic(:,:,:) ! phonon angular momentum per mode per atom

!! Ph+MBC Eigenvalue problem variables
COMPLEX(dp), ALLOCATABLE :: omega_q(:,:)

!! Lapack variables
COMPLEX(dp), ALLOCATABLE :: vl(:,:), vr(:,:)
COMPLEX(dp),ALLOCATABLE :: matA(:,:),matB(:,:)
REAL(dp), ALLOCATABLE :: w(:)
COMPLEX(dp), ALLOCATABLE :: w_cplx(:),alpha(:),beta(:)
INTEGER, ALLOCATABLE :: ipiv(:)
INTEGER :: n

!! Phangmom variables
REAL(dp) :: phangmom_net(3) ! net phonon angular momentum
REAL(dp),ALLOCATABLE :: phangmom_net_temp(:,:), temp(:)
REAL(dp) :: temp_min, temp_max, occ
INTEGER :: temp_n
REAL(dp),ALLOCATABLE :: atomic_character(:,:)
REAL(dp),ALLOCATABLE :: phmag_character(:,:)

!! Visualization variables
character(len=len100) :: prefix='outdata/dump_'
REAL(dp),ALLOCATABLE :: traj(:,:,:)
REAL(dp) :: timestep
INTEGER :: nperiod
! Variables END

call readInput(crystName)
if(trim(crystName)=='CrI3') then
  natom_species=2
  ALLOCATE(atomic_masses(natom_species))
  ALLOCATE(atom_species_label(natom_species))
  atomic_masses(:)=[51.9961_dp, 126.90447_dp]
  atom_species_label(:)=['Cr','I ']

  natom = 8
  ALLOCATE(atom_species(natom))
  atom_species(:)=[1,1,2,2,2,2,2,2]

  nmagatom=2
  ALLOCATE(magatom(nmagatom))
  magatom(:)=[1,2]
  nmagdir=2 ! only x and y directions

  ALLOCATE(xred(3,natom))
  xred(:,:)=reshape(\
      [  0.166280508349_dp,      0.166280508086_dp,      0.166280503996_dp,\
        -0.166280505848_dp,     -0.166280506097_dp,     -0.166280506002_dp,\
         0.424051945955_dp,     -0.220361170953_dp,      0.065302026100_dp,\
         0.065302024331_dp,      0.424051943912_dp,     -0.220361167348_dp,\
        -0.220361171173_dp,      0.065302034889_dp,      0.424051945241_dp,\
        -0.424051832996_dp,      0.220361112829_dp,     -0.065302080744_dp,\
        -0.065302086479_dp,     -0.424051832789_dp,      0.220361116534_dp,\
         0.220361117860_dp,     -0.065302089878_dp,     -0.424051837778_dp],\
      shape(xred))
  latt_vec(:,:)=reshape(\
      [ 7.4921305100422115_dp, 0.0_dp,               12.476606612745105_dp,\
       -3.7460652550211058_dp, 6.488375350135703_dp, 12.476606612745105_dp,\
       -3.7460652550211058_dp,-6.488375350135703_dp, 12.476606612745105_dp],\
       shape(latt_vec))
else if(trim(crystName)=='CrBr3') then
    natom_species=2
  ALLOCATE(atomic_masses(natom_species))
  ALLOCATE(atom_species_label(natom_species))
  atomic_masses(:)=[51.9961_dp, 79.904_dp]
  atom_species_label(:)=['Cr','Br']

  natom = 8
  ALLOCATE(atom_species(natom))
  atom_species(:)=[1,1,2,2,2,2,2,2]

  nmagatom=2
  ALLOCATE(magatom(nmagatom))
  magatom(:)=[1,2]
  nmagdir=2 ! only x and y directions

  ALLOCATE(xred(3,natom))
  xred(:,:)=reshape(\
      [  0.1662882979_dp,    0.1662882979_dp,    0.1662882979_dp,\
         0.8337117021_dp,    0.8337117021_dp,    0.8337117021_dp,\
         0.4216195246_dp,    0.7702368596_dp,    0.0694560195_dp,\
         0.0694560195_dp,    0.4216195246_dp,    0.7702368596_dp,\
         0.7702368596_dp,    0.0694560195_dp,    0.4216195246_dp,\
         0.5783804754_dp,    0.2297631404_dp,    0.9305439805_dp,\
         0.9305439805_dp,    0.5783804754_dp,    0.2297631404_dp,\
         0.2297631404_dp,    0.9305439805_dp,    0.5783804754_dp],\
      shape(xred))
  latt_vec(:,:)=reshape(\
      [  6.7493259557_dp,    0.0000000000_dp,   11.1229631832_dp,\
        -3.3746629780_dp,    5.8450877362_dp,   11.1229631832_dp,\
        -3.3746629780_dp,   -5.8450877362_dp,   11.1229631832_dp],\
       shape(latt_vec))
end if

nat3=3*natom
nat6=6*natom
nmagpert=nmagatom*nmagdir

WRITE(*,*) "Initialization of basic quantities"
ALLOCATE(xcart(3,natom))
ALLOCATE(amus(natom))
do iat=1,natom
  amus(iat) = atomic_masses(atom_species(iat))
  xcart(:,iat) = matmul(latt_vec,xred(:,iat))
end do


WRITE(*,*) "Read dynamical matrix from file (",fname_kpp,")"
WRITE(*,*) "Convention is: cartesian coordinates"
WRITE(*,*) "Units: Ha/Bohr^2" ! but read as eV/Angstrom^2
ALLOCATE(k_pp(3*natom,3*natom))
open(newunit=io, file=fname_kpp)
do irow=1,3*natom
  read(io,*) k_pp(irow,1:3*natom)
  k_pp(irow,1:3*natom) = k_pp(irow,1:3*natom) / (27.2114_dp/(0.52917_dp)**2)
end do
close(io)
WRITE(*,*) "Read"

WRITE(*,*) "Applying ASR on reciprocal IFCs"
call apply_asr_dynmat_q0(k_pp,natom)

WRITE(*,*) "Applying Hermitian symmetry to the IFCs"
call symmetrize_hermitian(k_pp,3*natom)

!Write the results
WRITE(*,*) ' Reciprocal interatomic force constants (Ha/Bohr^2)'
WRITE(*,*) '  atom1  dir  atom2  dir        Real              Imag'
ii=1
outer_loop:\
do iat1=1,natom
  do idir1=1,3
    ipc1=(iat1-1)*3 + idir1
    do iat2=1,natom
      do idir2= 1, 3
        ipc2=(iat2-1)*3 + idir2
        WRITE(*,'(2(i4,4x,a2,2x),2x,2es18.9)') &
        & iat1, cart(idir1), iat2, cart(idir2), &
        & real(k_pp(ipc1,ipc2)), aimag(k_pp(ipc1,ipc2))
        ii = ii + 1
        if(prt_debug/=1 .and. ii>max_prt_lines) then
          WRITE(*,*) "Enough output, stop here."
          EXIT outer_loop
        end if
      end do
    end do
    WRITE(*,*) '   '
  end do
end do outer_loop

WRITE(*,*) "Divide IFC by the square root of the masses"
call massmult(k_pp, natom, amus)

WRITE(*,*)
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,'(a)') "-------------------               Bare Phonons               -------------------"
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,*)

ALLOCATE(tmp_matrix(nat3,nat3))
ALLOCATE(w(nat3))
tmp_matrix(:,:) = k_pp(:,:) ! copy otherwise it'll be destroyed

WRITE(*,*) char(10) // " Solving the eigenvalue problem..."
call zheev_wrp(tmp_matrix, nat3, w)

! copy the eigenvectors before erasing them inadvertedly
ALLOCATE(pheigvec(3,natom,nat3))
pheigvec(:,:,:) = RESHAPE(tmp_matrix,[3,natom,nat3])
DEALLOCATE(tmp_matrix)

ALLOCATE(phfreq(nat3))
do ii=1,nat3
  if (w(ii) > tol12) then
    phfreq(ii) = SQRT(w(ii))
  else if (w(ii) > -tol12) then
    phfreq(ii) = zero
 else
    phfreq(ii) = -SQRT(-w(ii))
  end if
end do
DEALLOCATE(w)

WRITE(*,*) " Phonon frequencies (in cm-1 and in meV):"
WRITE(*,'(a)') "Mode, frequency(cm-1),   frequency(meV)"
do ii=1,3*natom
  WRITE(*,'(i3,2x,es15.8,4x,f15.10)') ii, phfreq(ii)*Ha_2_cm, phfreq(ii)*Ha_2_meV
end do

if(prt_debug==1) then
  WRITE(*,*) " Printing the displacement"
  do ii=1,nat3
    WRITE(*,'(a,i3)') "Mode: ", ii
    WRITE(*,'(a)') " Re(x)        Im(x)          Re(y)        Im(y)          Re(z)        Im(z)"
    do jj=1,natom
      WRITE(*,'(3(es11.4,2x,es11.4,4x))') pheigvec(1:3,jj,ii)
    end do
    WRITE(*,*) ""
  end do
end if

WRITE(*,*) " Examining the angular momentun of the phonon modes..."
WRITE(*,*) " Everything should be zero..."
ALLOCATE(phangmom(3,nat3))
phangmom(:,:) = 0
do imode=1,nat3
  do iat=1,natom
    phangmom(1,imode) = phangmom(1,imode) +\
                        2 * ( real(pheigvec(2,iat,imode)) * aimag(pheigvec(3,iat,imode))\
                            -aimag(pheigvec(2,iat,imode)) *  real(pheigvec(3,iat,imode)))
    phangmom(2,imode) = phangmom(2,imode) +\
                        2 * ( real(pheigvec(3,iat,imode)) * aimag(pheigvec(1,iat,imode))\
                            -aimag(pheigvec(3,iat,imode)) *  real(pheigvec(1,iat,imode)))
    phangmom(3,imode) = phangmom(3,imode) +\
                        2 * ( real(pheigvec(1,iat,imode)) * aimag(pheigvec(2,iat,imode))\
                            -aimag(pheigvec(1,iat,imode)) *  real(pheigvec(2,iat,imode)))
  end do
end do

WRITE(*,'(a)') "Phonon angular momentum (in units of hbar):"
do imode=1,nat3
  WRITE(*,'(1x,i2,2x,3(es18.9,2x))') imode, phangmom(:,imode)
end do

DEALLOCATE(pheigvec)
DEALLOCATE(phangmom)


WRITE(*,*)
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,'(a)') "-------------------         Phonons+mbc (Saparov22)          -------------------"
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,*)

WRITE(*,*) "We follow the resolution detailed in"
WRITE(*,*) "Saparov et al., Phys. Rev. B 105, no. 6 (2022): 064303."
WRITE(*,*) "https://doi.org/10.1103/PhysRevB.105.064303."
WRITE(*,*)

WRITE(*,*) "Read mbc matrix from file (",fname_gpp,")"
WRITE(*,*) "Convention is: cartesian coordinates"
WRITE(*,*) "Units: 1/Bohr^2" ! but read as 1/Angstrom^2
ALLOCATE(g_pp(3*natom,3*natom))
open(newunit=io, file=fname_gpp)
do irow=1,3*natom
  read(io,*) g_pp(irow,1:3*natom)
  g_pp(irow,1:3*natom) = g_pp(irow,1:3*natom) / (1/(0.52917_dp)**2)
end do
close(io)
WRITE(*,*) "Read"

! WRITE(*,*) "Applying ASR on reciprocal IVFCs"
! call apply_asr_dynmat_q0(g_pp,natom)

WRITE(*,*) "Applying anti-Hermitian symmetry to the MBC"
call symmetrize_antihermitian(g_pp,3*natom)

!Write the results
WRITE(*,*) ' Molecular Berry Curvature (1/Bohr^2)'
WRITE(*,*) '  atom1  dir  atom2  dir        Real              Imag'
ii=1
outer_loop2:\
do iat1=1,natom
  do idir1=1,3
    ipc1=(iat1-1)*3 + idir1
    do iat2=1,natom
      do idir2= 1, 3
        ipc2=(iat2-1)*3 + idir2
        WRITE(*,'(2(i4,4x,a2,2x),2x,2es18.9)') &
        & iat1, cart(idir1), iat2, cart(idir2), &
        & real(g_pp(ipc1,ipc2)), aimag(g_pp(ipc1,ipc2))
        ii = ii + 1
        if(prt_debug/=1 .and. ii>max_prt_lines) then
          WRITE(*,*) "Enough output, stop here."
          EXIT outer_loop2
        end if
      end do
    end do
    WRITE(*,*) '   '
  end do
end do outer_loop2

WRITE(*,*) "Divide MBC by the square root of the masses"
call massmult(g_pp, natom, amus)

! call solve_saparov(natom,k_pp,g_pp,phfreq,prtdebug)
! This was initially meant as a test, disabled for now
! You can reenable it, if you break things and want to have a reference implementation
WRITE(*,*) "We skip the Saparov resolution (uncomment line code if you really want it)"

WRITE(*,*)
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,'(a)') "-------------------             Phonons+mbc (Mignolet26)              ----------"
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,*)

! This was a first method, it doesn't match the improved formalism of Mignolet26 anymore
! mu is pre-multiplied by -i here
! However it gives the same results
! System to solve:
!  w (eps) = (1j*G_q,    I  ) (eps)
!    (mu )   (D_q,    1j*G_q) (mu )
! constructing big matrix:
ALLOCATE(omega_q(6*natom,6*natom))
omega_q(:,:) = zero
omega_q(1:nat3,1:nat3) = j_dpc*g_pp
do ii=1,nat3
  omega_q(ii,nat3+ii) = one
end do
omega_q(nat3+1:nat6,nat3+1:nat6) = j_dpc*g_pp
omega_q(nat3+1:nat6,1:nat3) = k_pp + MATMUL(TRANSPOSE(CONJG(g_pp)), g_pp)

ALLOCATE(w_cplx(nat6))
ALLOCATE(eigvec(nat6,nat6))
call zgeev_vr_wrp(omega_q, nat6, w_cplx, eigvec)

! TODO: Improve intelligence of this check
if (ANY(AIMAG(w_cplx) > tol5)) then
  WRITE(*,*) 'Imaginary freq detected, stopping here'
  STOP
end if
ALLOCATE(w(nat6))
w(:) = REAL(w_cplx(:))
DEALLOCATE(w_cplx)

! Sorting modes according to their frequencies
! as a word of caution the zero frequency modes are not necessarily correctly
! sorted. If you need to reuse them, be careful..
ALLOCATE(ind(nat6))
call sort_index(w, ind)

WRITE(*,*) "Newly computed phonon frequencies | correction (meV): "
do imode=1,nat3
  WRITE(*,*) imode, w(nat3+imode) * Ha_2_meV, (w(nat3+imode)-phfreq(imode)) * Ha_2_meV
end do

! Sorting the eigenvectors
! Using vl as a work array
ALLOCATE(vl(nat6,nat6))
do imode=1,nat6
  vl(:,imode) = eigvec(:,ind(imode))
end do
eigvec(:,:) = vl(:,:)
DEALLOCATE(vl)
DEALLOCATE(ind)
ALLOCATE(pheigvec(3,natom,nat3))
ALLOCATE(pheigvec_vel(3,natom,nat3))
pheigvec(:,:,:)     = RESHAPE(eigvec(1:nat3,     nat3+1:nat6),shape=[3,natom,nat3])
pheigvec_vel(:,:,:) = RESHAPE(eigvec(nat3+1:nat6,nat3+1:nat6),shape=[3,natom,nat3])
call renorm_phmbc(pheigvec, pheigvec_vel, natom, nmode=3*natom)

if(prt_debug==1) then
  do ii=1,nat3
    WRITE(*,'(a,i3)') "Mode: ", ii
    WRITE(*,'(a)') " Re(x)        Im(x)          Re(y)        Im(y)          Re(z)        Im(z)"
    do jj=1,natom
      WRITE(*,'(3(es11.4,2x,es11.4,4x))') pheigvec(1:3,jj,ii) ! /sqrt(amus(jj)*amu_emass)
    end do
    ! ! velocity part
    ! WRITE(*,'(a)') "Velocity part"
    ! do jj=1,natom
    !   WRITE(*,'(3(es11.4,2x,es11.4,4x))') pheigvec_vel(1:3,jj,ii)
    ! end do
    ! ! comparison
    ! WRITE(*,'(a)') "Displ / Velocity"
    ! do jj=1,natom
    !   WRITE(*,'(3(es11.4,2x,es11.4,4x))') pheigvec(1:3,jj,ii) / pheigvec_vel(1:3,jj,ii)
    ! end do
    WRITE(*,*) ""
  end do
end if

! Computing the angular momentum
ALLOCATE(phangmom_atomic(3,natom,nat3))
call compute_phangmom_phmbc(pheigvec, pheigvec_vel, phangmom_atomic, natom, nat3)
ALLOCATE(phangmom(3,nat3))
phangmom(:,:) = sum(phangmom_atomic(:,:,:),dim=2) ! sum over atomic contributions

WRITE(*,'(a)') "Phonon angular momentum (in units of hbar and in cart. coord):"
do imode=1,nat3
  WRITE(*,'(1x,i2,2x,3(es18.9,2x))') imode, phangmom(:,imode)
end do

! Computing the amount of angular momentum carried by phonons
! L_net = sum_\nu L_\nu (1/2 + f(omega_\nu))
! at 0K we have L_net = sum_\nu L_\nu * 1/2
! This only accounts for Gamma phonons though
phangmom_net(:) = zero
do imode=4,nat3
  phangmom_net(:) = phangmom_net(:) + phangmom(:,imode) * 1/2
end do

WRITE(*,*) ""
WRITE(*,'(a)') "Net phonon angular momentum at 0K (in units of hbar and in cart. coord): "
WRITE(*,'(3(es18.9,2x))') phangmom_net(:)

! Now let's look at finite temperature
! Let's go from 0 to 1000K by steps of 50K
temp_n=21
temp_min=0
temp_max=1000
ALLOCATE(temp(temp_n))
do ii=1,temp_n
  temp(ii) = temp_min + (ii-1)*(temp_max-temp_min)/(temp_n-1)
end do

ALLOCATE(phangmom_net_temp(3,temp_n))
phangmom_net_temp(:,:) = zero

! do the computations
do ii=1,temp_n
  do imode=4,nat3
    occ = occ_bose_einstein(w(nat3+imode), temp(ii))
    phangmom_net_temp(:,ii) = phangmom_net_temp(:,ii) + phangmom(:,imode) * (half + occ)
  end do
end do

! Output to screen
WRITE(*,*) ""
WRITE(*,'(a)') "Net phonon angular momentum as a function of temperature (in units of hbar and",\
               "in cart. coord): "
WRITE(*,'(a)') "Remark: the acoustic modes are not taken into account"

do ii=1,temp_n
  WRITE(*,'(f7.1,3x,3(es18.9,2x))') temp(ii), phangmom_net_temp(:,ii)
end do

ALLOCATE(atomic_character(natom_species,nat3))
call compute_atomic_charac_per_phmbc_mode(pheigvec,pheigvec_vel,natom,nat3,\
                                atom_species,natom_species,atomic_character)

WRITE(*,*)
WRITE(*,*) "Compute atomic contribution to each mode"
WRITE(*,'(a)') "imode     atom_type1        atom_type2        Total"
do imode=1,nat3
  WRITE(*,'(i3,2x,f16.9,2x,f16.9,4x,f14.9)') imode, atomic_character(1:2,imode), sum(atomic_character(:,imode))
end do

! set to true if you want some visualization
if (.False.) then
  ! Let us generate the trajectories of the atoms for visualization
  WRITE(*,'(a)') "Visualization..."
  do imode=15,16
    WRITE(*,'(a,i4)') ' Processing mode: ', imode
    timestep=1d-1 ! in picosec -> 100 femtosec
    nperiod=5
    call gen_ph_trajectory(pheigvec(:,:,imode),phfreq(imode),xcart,natom,amus,\
                           timestep,nperiod,traj)

    WRITE(msg,'(a,i2,a)') trim(prefix),imode,'_'
    call dump_ph_trajectory(msg,latt_vec,natom,atom_species,natom_species,atom_species_label,traj)
  end do
  WRITE(*,'(a)') "Done"
end if

DEALLOCATE(pheigvec)
DEALLOCATE(phangmom_atomic)
DEALLOCATE(phangmom)
DEALLOCATE(w)


WRITE(*,*)
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,'(a)') "-------------------           Bare Magnons (Niu98)           -------------------"
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,*)

WRITE(*,*) "We follow the resolution detailed in"
WRITE(*,*) "Niu and Kleinman, Phys. Rev. Lett. 80, no. 10 (1998): 2205–8."
WRITE(*,*) "https://doi.org/10.1103/PhysRevLett.80.2205."
WRITE(*,*)

WRITE(*,*) "Read spin-spin K_ss matrix from file (",fname_kss,")"
WRITE(*,*) "Convention is: cartesian coordinates"
WRITE(*,*) "Units: Ha" ! but read as eV
ALLOCATE(k_ss(nmagpert,nmagpert))
open(newunit=io, file=fname_kss)
do irow=1,nmagpert
  read(io,*) k_ss(irow,1:nmagpert)
  k_ss(irow,1:nmagpert) = k_ss(irow,1:nmagpert) / 27.2114_dp
end do
close(io)
WRITE(*,*) "Read"

WRITE(*,*) "Read spin-spin Berry Curv. matrix from file (",fname_gss,")"
WRITE(*,*) "Convention is: cartesian coordinates"
WRITE(*,*) "Units: -"
ALLOCATE(g_ss(nmagpert,nmagpert))
open(newunit=io, file=fname_gss)
do irow=1,nmagpert
  read(io,*) g_ss(irow,1:nmagpert)
end do
close(io)
WRITE(*,*) "Read"

! Maybe apply ASR and other symm.
call symmetrize_hermitian(k_ss, nmagpert)
call symmetrize_antihermitian(g_ss, nmagpert)

!Write the results
WRITE(*,*) ' Inverse local spin susceptibility (Ha)'
WRITE(*,*) '  atom1  dir  atom2  dir        Real              Imag'
ii=1
ipc1=1
outer_loop3:\
do iat1=1,nmagatom
  do idir1=1,nmagdir
    ipc2=1
    do iat2=1,nmagatom
      do idir2= 1,nmagdir
        WRITE(*,'(2(i4,4x,a2,2x),2x,2es18.9)') &
        & magatom(iat1), cart(idir1), magatom(iat2), cart(idir2), &
        & real(k_ss(ipc1,ipc2)), aimag(k_ss(ipc1,ipc2))
        ii = ii + 1
        if(prt_debug/=1 .and. ii>max_prt_lines) then
          WRITE(*,*) "Enough output, stop here."
          EXIT outer_loop3
        end if
        ipc2=ipc2+1
      end do
    end do
    WRITE(*,*) '   '
    ipc1=ipc1+1
  end do
end do outer_loop3

!Write the results
WRITE(*,*) ' Spin Berry Curvature (-)'
WRITE(*,*) '  atom1  dir  atom2  dir        Real              Imag'
ii=1
ipc1=1
outer_loop4:\
do iat1=1,nmagatom
  do idir1=1,nmagdir
    ipc2=1
    do iat2=1,nmagatom
      do idir2=1,nmagdir
        WRITE(*,'(2(i4,4x,a2,2x),2x,2es18.9)') &
        & magatom(iat1), cart(idir1), magatom(iat2), cart(idir2), &
        & real(g_ss(ipc1,ipc2)), aimag(g_ss(ipc1,ipc2))
        ii = ii + 1
        if(prt_debug/=1 .and. ii>max_prt_lines) then
          WRITE(*,*) "Enough output, stop here."
          EXIT outer_loop4
        end if
        ipc2=ipc2+1
      end do
    end do
    WRITE(*,*) '   '
    ipc1=ipc1+1
  end do
end do outer_loop4

! Solving for the bare magnons: +/- iw G s = K s
WRITE(*,'(a)') "Solving system..."
ALLOCATE(matA(nmagpert,nmagpert))
ALLOCATE(matB(nmagpert,nmagpert))
ALLOCATE(alpha(nmagpert))
ALLOCATE(beta(nmagpert))
ALLOCATE(vr(nmagpert,nmagpert))
matA(:,:) = k_ss
matB(:,:) = g_ss
call zggev3_vr_wrp(matA, matB, nmagpert, alpha, beta, vr)
DEALLOCATE(matA)
DEALLOCATE(matB)
WRITE(*,'(a)') "Solved."

WRITE(*,'(a)') "Bare magnon frequencies (meV):"
do i=1,nmagpert
  if (abs(beta(i)) < tol12) then
    WRITE(*,*) "alhpa, beta: ", alpha(i), beta(i)
  end if
  WRITE(*,'(i3,2x,f16.10,a,f16.10)') i, AIMAG(alpha(i)/beta(i))*Ha_2_meV, ' Real part: ', REAL(alpha(i)/beta(i))*Ha_2_meV
end do

do ii=1,nmagpert
  WRITE(*,'(a,i3)') "Mode: ", ii
  WRITE(*,'(a)') " Re(s_x)      Im(s_x)        Re(s_y)      Im(s_y)"
  do jj=1,nmagatom
    WRITE(*,'(2(es11.4,2x,es11.4,4x))') real(vr(nmagdir*(jj-1)+1,ii)), aimag(vr(nmagdir*(jj-1)+1,ii)),&
                                        real(vr(nmagdir*(jj-1)+2,ii)), aimag(vr(nmagdir*(jj-1)+2,ii))
  end do
  WRITE(*,*) ""
end do
DEALLOCATE(vr)

WRITE(*,*)
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,'(a)') "-------------------          Phonons+Magnons (Ren24)         -------------------"
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,*)
! Solving method from Ren S., PRX 14 (1) 011041 (2024)

WRITE(*,*) "We follow the resolution detailed in"
WRITE(*,*) "Ren et al., Phys. Rev. X 14, no. 1 (2024): 011041."
WRITE(*,*) "https://doi.org/10.1103/PhysRevX.14.011041."
WRITE(*,*)

! Reading all the files
WRITE(*,*) "Read files"
WRITE(*,*) "Reading k_pp..." ! Units: Ha/Bohr^2, but read as eV/Angstrom^2
! ALLOCATE(k_pp(3*natom,3*natom))
open(newunit=io, file=fname_kpp)
do irow=1,3*natom
  read(io,*) k_pp(irow,1:3*natom)
  k_pp(irow,1:3*natom) = k_pp(irow,1:3*natom) / (27.2114_dp/(0.52917_dp)**2)
end do
call apply_asr_dynmat_q0(k_pp,natom)
call symmetrize_hermitian(k_pp,3*natom)
close(io)

WRITE(*,*) "Reading k_ps..." ! Units: Ha/Bohr, but read as eV/Angstrom
ALLOCATE(k_ps(3*natom,nmagpert))
open(newunit=io, file=fname_kps)
do irow=1,3*natom
  read(io,*) k_ps(irow,1:nmagpert)
  k_ps(irow,1:nmagpert) = k_ps(irow,1:nmagpert) / (27.2114_dp/0.52917_dp)
end do
close(io)
ALLOCATE(k_sp(nmagpert,nat3))
k_sp(:,:) = TRANSPOSE(CONJG(k_ps))

WRITE(*,*) "Reading k_ss..." ! Units: Ha, but read as eV
! ALLOCATE(k_ss(nmagpert,nmagpert))
open(newunit=io, file=fname_kss)
do irow=1,nmagpert
  read(io,*) k_ss(irow,1:nmagpert)
  k_ss(irow,1:nmagpert) = k_ss(irow,1:nmagpert) / 27.2114_dp
end do
call symmetrize_hermitian(k_ss, nmagpert)
close(io)

WRITE(*,*) "Reading g_pp..." ! Units: 1/Bohr^2, but read as 1/Angstrom^2
! ALLOCATE(g_pp(3*natom,3*natom))
open(newunit=io, file=fname_gpp)
do irow=1,3*natom
  read(io,*) g_pp(irow,1:3*natom)
  g_pp(irow,1:3*natom) = g_pp(irow,1:3*natom) / (1/(0.52917_dp)**2)
end do
call symmetrize_antihermitian(g_pp,3*natom)
close(io)

WRITE(*,*) "Reading g_ps..." ! Units: 1/Bohr, but read as 1/Angstrom
ALLOCATE(g_ps(3*natom,nmagpert))
open(newunit=io, file=fname_gps)
do irow=1,3*natom
  read(io,*) g_ps(irow,1:nmagpert)
  g_ps(irow,1:nmagpert) = g_ps(irow,1:nmagpert) / (1/0.52917_dp)
end do
close(io)
ALLOCATE(g_sp(nmagpert,nat3))
g_sp(:,:) = - TRANSPOSE(CONJG(g_ps))

WRITE(*,*) "Reading g_ss..." ! Units: -
open(newunit=io, file=fname_gss)
do irow=1,nmagpert
  read(io,*) g_ss(irow,1:nmagpert)
end do
call symmetrize_antihermitian(g_ss, nmagpert)
close(io)

! Solving method from Ren S., PRX 14 (1) 011041 (2024). https://doi.org/10.1103/PhysRevX.14.011041.
! This is the system to be solved
! dot{u}  = (0    1    0) (u)
! ddot{u} = (A_1 A_2 A_3) (dot{u})
! dot{s}  = (B_1 B_2 B_3) (s)
! With:
! A_1 = -M^-1 K_pp    + M^-1 G_ps G_ss^-1 K_sp
! A_2 = -M^-1 G_pp    - M^-1 G_ps G_ss^-1 G_sp
! A_3 = -M^-1 K_ps    + M^-1 G_ps G_ss^-1 K_ss
! B_1 =  G_ss^-1 K_sp
! B_2 = -G_ss^-1 G_sp
! B_3 =  G_ss^-1 K_ss
! Total size of the big matrix: 6*natom + nmagpert
! We need M^-1 and G_ss^-1 first. M^-1 is just (M^-1)_ij = diag(1/M_i)
! G_ss has to be inverted through zgetrf and zgetri.
! Afterwards we can precompute M^-1 G_ps G_ss^-1. And we finish the construction
! Of the main matrix on the fly and a call to zgeev to obtain the eigenvals/vecs
! Afterwards we will talk about normalization, phangmom and so on...

! Inversion of G_ss^-1
ALLOCATE(ipiv(nmagpert))
call zgetrf_wrp(g_ss, nmagpert, ipiv)
call zgetri_wrp(g_ss, nmagpert, ipiv)
DEALLOCATE(ipiv)

! Precomputation of M^-1 G_ps G_ss^-1
! size of the matrix: nat3 x nmagpert
ALLOCATE(tmp_matrix(nat3,nmagpert))
call massmult_diag(g_ps, natom, amus) ! G_ps = M^-1 G_ps
tmp_matrix(:,:) = MATMUL(g_ps, g_ss) ! tmp_mat = M^-1 G_ps G_ss^-1

! Construction of the big matrix
WRITE(*,*) "Construction of the big matrix..."
n = nat6 + nmagpert
ALLOCATE(matA(n,n))
matA(:,:) = zero
! First row
matA(1:nat3,:) = zero
do i=1,nat3
  matA(i,nat3+i) = cone
end do
! Forming A_1
call massmult_diag(k_pp, natom, amus) ! K_pp = M^-1 K_pp
matA(nat3+1:nat6,1:nat3) = - k_pp + MATMUL(tmp_matrix, k_sp)
! Forming A_2
call massmult_diag(g_pp, natom, amus) ! G_pp = M^-1 G_pp
matA(nat3+1:nat6,nat3+1:nat6) = - 2*g_pp - MATMUL(tmp_matrix, g_sp)
! Forming A_3
call massmult_diag(k_ps, natom, amus) ! K_ps = M^-1 K_ps
matA(nat3+1:nat6,nat6+1:nat6+nmagpert) = - k_ps + MATMUL(tmp_matrix, k_ss)
! Forming B_1
matA(nat6+1:nat6+nmagpert,1:nat3) = MATMUL(g_ss, k_sp)
! Forming B_2
matA(nat6+1:nat6+nmagpert,nat3+1:nat6) = - MATMUL(g_ss, g_sp)
! Forming B_3
matA(nat6+1:nat6+nmagpert,nat6+1:nat6+nmagpert) = MATMUL(g_ss, k_ss)

! Find eigenvals/vecs
WRITE(*,*) "Solving the eigen problem..."
ALLOCATE(w_cplx(n))
ALLOCATE(vr(n,n))
call zgeev_vr_wrp(matA, n, w_cplx, vr)
DEALLOCATE(tmp_matrix)
DEALLOCATE(matA)
WRITE(*,*) "Solved."

! Retrieving imaginary part which corresponds to the modes frequency
ALLOCATE(w(n))
w(:) = AIMAG(w_cplx(:))
DEALLOCATE(w_cplx)

! Sorting modes according to their frequencies
ALLOCATE(ind(n))
call sort_index(w, ind)

! Sorting the eigenvectors
! Using vl as a work array
ALLOCATE(vl(n,n))
do imode=1,n
  vl(:,imode) = vr(:,ind(imode))
end do
vr(:,:) = vl(:,:)
DEALLOCATE(vl)
DEALLOCATE(ind)

WRITE(*,*) "Eigenvalues (meV):"
do ii=1,n/2 ! Only positive freq modes
  imode=n/2+ii
  WRITE(*,*) ii, w(imode) * Ha_2_meV
end do

! renormalize eigenvectors
call renorm_phmag_Ren24(vr,natom,n)

! Compute magnon-phonon character
ALLOCATE(phmag_character(2,n))
call compute_phmag_charac_per_mode_SRen24(vr, natom, nmagpert, n, phmag_character)
WRITE(*,*) "Phononic character:                Magnonic character:"
do ii=1,n/2 ! Only positive freq modes
  imode=n/2+ii
  WRITE(*,*) ii, phmag_character(:,imode)
end do
WRITE(*,*) "Total:"
WRITE(*,*) SUM(phmag_character, DIM=2)
WRITE(*,*) "Expected: ", nat6, nmagpert

DEALLOCATE(phmag_character)

WRITE(*,*)
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,'(a)') "-------------------       Phonons+Magnons (Mignolet26)       -------------------"
WRITE(*,'(a)') "--------------------------------------------------------------------------------"
WRITE(*,*)

WRITE(*,*) "Implemented but some cleaning should be done"

call solve_phmag(natom,nmagpert,amus,natom_species,atom_species,&
                &fname_kpp,fname_kss,fname_kps,&
                &fname_gpp,fname_gss,fname_gps)

WRITE(*,*) "Finalize..."
! Final frees
WRITE(*,*) "Deallocating common arrays."
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
DEALLOCATE(omega_q)
DEALLOCATE(phfreq)
DEALLOCATE(amus)

WRITE(*,*) 'Run for ', trim(crystName), ' finished.'

end program
