MODULE m_dynmat

! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

use m_commons

implicit none


private

public :: symmetrize_hermitian
public :: symmetrize_antihermitian
public :: massmult
public :: massmult_left
public :: massmult_right
public :: massmult_diag
public :: apply_asr_dynmat_q0

CONTAINS


subroutine symmetrize_hermitian(mat, npert)
!! FUNCTION
!! Symmetrize a Hermitian matrix
!!
!! INPUTS
!! mat = npertxnpert matrix to be Hermitianized
!!
!! OUTPUT
!! symMat = npertxnpert Hermitianized matrix

!Arguments ------------------------------------
integer, intent(in) :: npert
complex(dp),intent(inout) :: mat(npert,npert)
!Local variables ------------------------------
complex(dp),allocatable :: tmp(:,:)
! *************************************************************************

ALLOCATE(tmp(npert,npert))
tmp(:,:) = (mat + CONJG(TRANSPOSE(mat)))/2
mat(:,:) = tmp
DEALLOCATE(tmp)

end subroutine symmetrize_hermitian
!!***

subroutine symmetrize_antihermitian(mat, npert)
!! FUNCTION
!! Symmetrize an anti-Hermitian matrix
!!
!! INPUTS
!! mat = npertxnpert matrix to be Hermitianized
!!
!! OUTPUT
!! symMat = npertxnpert Hermitianized matrix

!Arguments ------------------------------------
integer, intent(in) :: npert
complex(dp),intent(inout) :: mat(npert,npert)
!Local variables ------------------------------
complex(dp),allocatable :: tmp(:,:)
! *************************************************************************

ALLOCATE(tmp(npert,npert))
tmp(:,:) = (mat - CONJG(TRANSPOSE(mat)))/2
mat(:,:) = tmp
DEALLOCATE(tmp)

end subroutine symmetrize_antihermitian
!!***

subroutine massmult(dmatrix, natom, amus)
!! FUNCTION
!! Multiplies a dynamical matrix by 1/sqrt(M_i*M_j)
!!
!! INPUTS
!! dmatrix(3*natom,3*natom): cplx*16
!! natom: number of atoms
!! amus(natom): atomic masses
!!
!! OUTPUT
!! dmatrix multiplied by 1/sqrt(M_i*M_j)

!Arguments ------------------------------------
integer,intent(in) :: natom
complex(dp),intent(inout) :: dmatrix(3*natom,3*natom)
real(dp),intent(in) :: amus(natom)
!Local variables ------------------------------
integer :: iat1,iat2,i,j,ii,jj
real(dp) :: m1,m2
! *************************************************************************

jj = 1
do iat2=1,natom
  m2 = amus(iat2)
  do j=1,3
    ii = 1
    do iat1=1,natom
      m1 = amus(iat1)
      do i=1,3
        dmatrix(ii,jj) = dmatrix(ii,jj) / SQRT(m1*m2) / amu_emass
        ii = ii + 1
      end do
    end do
    jj = jj +1
  end do
end do

end subroutine massmult
!!***

subroutine massmult_left(dmatrix, natom, dimB, amus)
!! FUNCTION
!! Multiplies a dynamical matrix by 1/sqrt(M_i)
!!
!! INPUTS
!! dmatrix(3*natom,dimB): cplx*16
!! natom: number of atoms
!! dimB: second dimension of the matrix
!! amus(natom): atomic masses
!!
!! OUTPUT
!! dmatrix multiplied by 1/sqrt(M_i)

!Arguments ------------------------------------
integer,intent(in) :: natom,dimB
complex(dp),intent(inout) :: dmatrix(3*natom,dimB)
real(dp),intent(in) :: amus(natom)
!Local variables ------------------------------
integer :: iat1,i,ii,jj
! *************************************************************************

do jj=1,dimB
  ii = 1
  do iat1=1,natom
    do i=1,3
      dmatrix(ii,jj) = dmatrix(ii,jj) / SQRT(amus(iat1)*amu_emass)
      ii = ii + 1
    end do
  end do
end do

end subroutine massmult_left
!!***

subroutine massmult_right(dmatrix, dimA, natom, amus)
!! FUNCTION
!! Multiplies a dynamical matrix by 1/sqrt(M_j)
!!
!! INPUTS
!! dmatrix(dimA,3*natom): cplx*16
!! dimA: first dimension of the matrix
!! natom: number of atoms
!! amus(natom): atomic masses
!!
!! OUTPUT
!! dmatrix multiplied by 1/sqrt(M_j)

!Arguments ------------------------------------
integer,intent(in) :: dimA,natom
complex(dp),intent(inout) :: dmatrix(dimA,3*natom)
real(dp),intent(in) :: amus(natom)
!Local variables ------------------------------
integer :: ii,iat2,j,jj
! *************************************************************************

jj = 1
do iat2=1,natom
  do j=1,3
    do ii=1,dimA
      dmatrix(ii,jj) = dmatrix(ii,jj) / SQRT(amus(iat2)*amu_emass)
    end do
    jj = jj +1
  end do
end do

end subroutine massmult_right
!!***

subroutine massmult_diag(dmatrix, natom, amus)
!! FUNCTION
!! Multiplies on the left a dynamical matrix by delta_ij/M_i
!! B(i,j) = sum_k (M^-1)_ik * A(k,j)
!!        = 1/M_i * A(i,j)
!!
!! INPUTS
!! dmatrix(3*natom,:): cplx*16
!! natom: number of atoms
!! amus(natom): atomic masses
!!
!! OUTPUT
!! dmatrix multiplied by delta_ij/M_i

!Arguments ------------------------------------
integer,intent(in) :: natom
complex(dp),intent(inout) :: dmatrix(:,:)
real(dp),intent(in) :: amus(natom)
!Local variables ------------------------------
integer :: iat
! *************************************************************************

do iat=1,natom
  dmatrix(3*(iat-1)+1,:) = dmatrix(3*(iat-1)+1,:) / (amus(iat)*amu_emass)
  dmatrix(3*(iat-1)+2,:) = dmatrix(3*(iat-1)+2,:) / (amus(iat)*amu_emass)
  dmatrix(3*(iat-1)+3,:) = dmatrix(3*(iat-1)+3,:) / (amus(iat)*amu_emass)
end do

end subroutine massmult_diag
!!***

subroutine apply_asr_dynmat_q0(dynmat,natom)
!! FUNCTION
!! Apply the acoustic sum rule on the dynamical matrix at q=0
!!
!! INPUTS
!! dynmat(3,natom,3,natom): dynamical matrix at q=0, before dividing by the
!!  square roots of the masses.
!! natom: number of atoms
!!
!! OUTPUT
!! dynmat with asr applied
!! atmfrc(:,iat1,:,iat1) -= sum_{iat2} atmfrc(:,iat2,:,iat1)

!Arguments ------------------------------------
integer,intent(in) :: natom
complex(dp),intent(inout) :: dynmat(3,natom,3,natom)
!Local variables ------------------------------
integer :: iat1,iat2
complex(dp) :: corr(3,3,natom)
! *************************************************************************

corr = zero
do iat1=1,natom
  do iat2=1,natom
    corr(:,:,iat1) = corr(:,:,iat1) + dynmat(:,iat2,:,iat1)
  end do
  dynmat(:,iat1,:,iat1) = dynmat(:,iat1,:,iat1) - corr(:,:,iat1)
end do

end subroutine apply_asr_dynmat_q0
!!***

END MODULE
