! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

module m_solve_saparov

USE m_commons
USE m_lapack_wrappper


IMPLICIT NONE

PRIVATE

PUBLIC :: solve_saparov

contains


subroutine solve_saparov(natom,k_pp,g_pp,phfreq_bare,prt_debug)
!! FUNCTION
!! Compute phonon+mbc frequencies using Saparov's method:
!! Saparov D. et al., PRB 105 (6), 064303 (2022)
!!
!! INPUTS
!! k_pp: dynamical matrix 3natx3nat
!! g_pp: mbc        matrix 3natx3nat
!!
!! OUTPUT
!! Nothing
!!
!! SIDE EFFECT
!! Compute the phonon frequencies including the mbc

!Arguments ------------------------------------
integer,intent(in) :: natom
complex(dp),intent(in) :: k_pp(3*natom,3*natom)
complex(dp),intent(in) :: g_pp(3*natom,3*natom)
REAL(dp),intent(in) :: phfreq_bare(3*natom) ! bare phonon frequencies ! used only for comparison
integer,intent(in) :: prt_debug
!Local variables ------------------------------
! Scalar
INTEGER :: i,ii,jj,kk,iat,imode
INTEGER :: nat3, nat6
real(dp),PARAMETER :: omega_0 = 1.0D0/Ha_2_meV
REAL(dp) :: norm
! Array
COMPLEX(dp), ALLOCATABLE :: D_q(:,:)
COMPLEX(dp), ALLOCATABLE :: sigma_x(:,:)
COMPLEX(dp), ALLOCATABLE :: omega_q(:,:), omega_q_sqrt(:,:)
COMPLEX(dp), ALLOCATABLE :: H_q(:,:)
COMPLEX(dp),ALLOCATABLE :: eigvec(:,:) ! general eignvec
COMPLEX(dp),ALLOCATABLE :: pheigvec(:,:,:) ! phonon eigenvec(idir,iat,imode)
REAL(dp),ALLOCATABLE :: phangmom(:,:) ! phonon angular momentum per mode
! Lapack Arrays
REAL(dp), ALLOCATABLE :: w(:)
! Other Arrays

nat3=3*natom
nat6=6*natom

WRITE(*,*) " Constructing intermediate matrices..."
! construct sigma_x
ALLOCATE(sigma_x(6*natom,6*natom))
sigma_x = zero
do ii=1,3*natom
  sigma_x(3*natom+ii,ii) = one
  sigma_x(ii,3*natom+ii) = one
end do

! construct D_q = K_q + G_q^dagger * G_q
ALLOCATE(D_q(3*natom,3*natom))
D_q(:,:) = k_pp + MATMUL(TRANSPOSE(CONJG(g_pp)), g_pp)
! construct Omega_q
! Omega_q = (D_q,    -1j*G_q^T)
!           (1J*G_q,     I    )
ALLOCATE(omega_q(6*natom,6*natom))
omega_q(1:3*natom, 1:3*natom) = D_q / omega_0
omega_q(1:3*natom, 3*natom+1:6*natom) = -j_dpc*TRANSPOSE(CONJG(g_pp))
omega_q(3*natom+1:6*natom, 1:3*natom) = j_dpc*g_pp
do ii=3*natom+1,6*natom
  omega_q(ii,ii) = one * omega_0
end do
DEALLOCATE(D_q)

! construct Omega_q^1/2
! First we need the spectral decomposition of Omega_q: U Lambda U^-1
! Then Omega_q^1/2 is given by U Lambda^1/2 U^-1
ALLOCATE(w(nat6))
call zheevd_wrp(omega_q, nat6, w)
WRITE(*,*) w, size(w)
do kk=1,nat6
  if (w(kk) > tol8) then
    cycle
  else if (w(kk) > -tol8) then
    w(kk) = zero
  else
    WRITE(*,*) "There was some problem somewhere"
    WRITE(*,*) "I'm returning here"
    RETURN ! STOP
  end if
end do
WRITE(*,*) "Eigenvalues of Omega_q for check (meV):"
do kk=1,nat6
  WRITE(*,*) kk, w(kk) * Ha_2_meV
end do

ALLOCATE(omega_q_sqrt(nat6,nat6))
omega_q_sqrt = zero
do jj=1,nat6
  do ii=1,nat6
    do kk=1,nat6
      ! Omega^1/2 (i,j) = sum_k  eigvec_i(k) * eigval^1/2 * (eigvec_j(k))^*
      omega_q_sqrt(ii,jj) = omega_q_sqrt(ii,jj) + omega_q(ii,kk) * sqrt(w(kk)) * CONJG(omega_q(jj,kk))
    end do
  end do
end do
DEALLOCATE(omega_q)

WRITE(*,*) " Constructing the big matrix..."
! construct H_q to solve
! H_q = Omega_q_sart * sigma_x * Omega_q_sqrt
ALLOCATE(H_q(nat6,nat6))
H_q(:,:) = MATMUL(MATMUL(omega_q_sqrt, sigma_x), omega_q_sqrt)

WRITE(*,*) " Solving the system..."
! For reference the system we are solving here is from:
! Saparov D. et al., PRB 105 (6), 064303 (2022)
call zheev_wrp(H_q, nat6, w)

WRITE(*,*) " System solved!"
! Saving the eigenvectors before doing anything else
! Normally they are in the form: Psi^tilde = Omega^1/2 * Psi,
! with Psi(1:nat3)      = phonon eigendispl
! and  Psi(nat3+1:nat6) = phonon eigenvelocities
! Let's try to retrieve the Psi
ALLOCATE(eigvec(nat6,nat6))
eigvec(:,:) = H_q(:,:)
DEALLOCATE(H_q)
call zhesv_wrp(omega_q_sqrt, nat6, eigvec)
DEALLOCATE(omega_q_sqrt)
! eigvec now contains the Psi
! Renormalization of the eigenvectors
do ii=1,nat6
  norm = REAL(DOT_PRODUCT(eigvec(:,ii), MATMUL(sigma_x, eigvec(:,ii))),dp)
  eigvec(:,ii) = eigvec(:,ii) / norm
end do
DEALLOCATE(sigma_x)

WRITE(*,*) "Newly computed phonon frequencies | correction (meV): "
i = 1
do ii=nat3+1,nat6
  WRITE(*,*) i, w(ii) * Ha_2_meV, (w(ii)-phfreq_bare(i)) * Ha_2_meV
  i = i + 1
end do

! Examining the angular momentum of the phonon modes now
WRITE(*,*) " Examining the angular momentun of the phonon modes..."

! First we need to extract the phonon part from the full eigenvectors
ALLOCATE(pheigvec(3,natom,nat3))
pheigvec(:,:,:) = RESHAPE(eigvec(nat3+1:nat6,1:nat3),[3,natom,nat3])
! we select positive freq and the phonon part

! Let's print the zero freq mode for check
if(prt_debug==1) then
  do ii=1,nat3
    WRITE(*,'(a,i3)') "Mode: ", ii
    WRITE(*,'(a)') " Re(x)        Im(x)          Re(y)        Im(y)          Re(z)        Im(z)"
    do jj=1,natom
      WRITE(*,'(3(es11.4,2x,es11.4,4x))') pheigvec(1:3,jj,ii)*sqrt(2*w(nat3+ii)/omega_0)
    end do
    WRITE(*,*) ""
  end do
end if

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

DEALLOCATE(w)
DEALLOCATE(eigvec)
DEALLOCATE(pheigvec)
DEALLOCATE(phangmom)


end subroutine solve_saparov

end module m_solve_saparov
