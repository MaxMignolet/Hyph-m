MODULE m_lapack_wrappper

! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

use m_commons

implicit none


private

public :: zheev_wrp
public :: zheevd_wrp
public :: zhesv_wrp
public :: zgeev_vr_wrp
public :: zggev3_vr_wrp
public :: zgetrf_wrp
public :: zgetri_wrp

! Lapack functions
EXTERNAL zheev
EXTERNAL zheevd
EXTERNAL zhegv
EXTERNAL zhesv
EXTERNAL zgeev
EXTERNAL zggev3
EXTERNAL zgetrf
EXTERNAL zgetri

CONTAINS


subroutine err_handler(info)
!! FUNCTION
!! Handles Lapack return info value
!!
!! INPUTS
!! info
!!
!! OUTPUT
!! Stop execution if info/=0

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: info
!Local variables ------------------------------
! *************************************************************************

if (info/=0) then
  WRITE(*,'(a)') "Error during call to lapack zheev routine."
  WRITE(*,'(a,1x,i2)') "Error code: ", info
  WRITE(*,'(a)') "STOPPING"
  STOP
end if

end subroutine err_handler
!!***

subroutine zheev_wrp(matA, n, w)
!! FUNCTION
!! Call to zheev. Computes the eigenvalues and eigenvector of the cplx*16
!! square Hermitian matrix of size n. The matrix is destroyed during the call.
!!
!! INPUTS
!! matA: input matrix (cplx*16, Hemirtian)
!! n: size of the matrix
!!
!! OUTPUT
!! matA contains the eigenvectors
!! w contains the eigenvalues

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(INOUT), DIMENSION(n,n) :: matA
REAL(dp), INTENT(OUT), DIMENSION(n) :: w
!Local variables ------------------------------
CHARACTER(1) :: jobz, uplo
INTEGER :: lwork, lrwork, info
COMPLEX(dp), ALLOCATABLE :: work(:), rwork(:)
! *************************************************************************

! specification of the problem:
jobz = 'V' ! eigenvalues AND eigenvectors to be returned
uplo = 'U' ! consider the upper diagonal
lwork = max(1,2*n-1)
lrwork = max(1, 3*n-2)

ALLOCATE(rwork(lrwork))

! Query optimal work size
lwork = -1
ALLOCATE(work(1))
call zheev(jobz, uplo, n, matA, n, w, work, lwork, rwork, info)
call err_handler(info)
lwork = INT(work(1))

! Actual call to zheev
DEALLOCATE(work)
ALLOCATE(work(lwork))
call zheev(jobz, uplo, n, matA, n, w, work, lwork, rwork, info)
call err_handler(info)

DEALLOCATE(work)
DEALLOCATE(rwork)

end subroutine zheev_wrp
!!***

subroutine zheevd_wrp(matA, n, w)
!! FUNCTION
!! Call to zheevd. Computes the eigenvalues and eigenvector of the cplx*16
!! square Hermitian matrix of size n. It uses a divide and conquer algorithm.
!! The matrix is destroyed during the call.
!!
!! INPUTS
!! matA: input matrix (cplx*16, Hemirtian)
!! n: size of the matrix
!!
!! OUTPUT
!! matA contains the eigenvectors
!! w contains the eigenvalues

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(INOUT), DIMENSION(n,n) :: matA
REAL(dp), INTENT(OUT), DIMENSION(n) :: w
!Local variables ------------------------------
CHARACTER(1) :: jobz, uplo
INTEGER :: lwork, lrwork, liwork, info
COMPLEX(dp), ALLOCATABLE :: work(:), rwork(:), iwork(:)
! *************************************************************************

! specification of the problem:
jobz = 'V' ! eigenvalues AND eigenvectors to be returned
uplo = 'U' ! consider the upper diagonal

! Query optimal work size
lwork = -1
ALLOCATE( work(1))
ALLOCATE(rwork(1))
ALLOCATE(iwork(1))
call zheevd(jobz, uplo, n, matA, n, w, work, lwork, rwork, lrwork,&
            iwork, liwork, info)
call err_handler(info)
lwork  = INT( work(1))
lrwork = INT(rwork(1))
! liwork = INT(iwork(1)) ! for some reason the workspace query does not return a
                         ! correct value for liwork. Putting it manually...
liwork = 3 + 5*n

! Actual call to zheev
DEALLOCATE( work)
DEALLOCATE(rwork)
DEALLOCATE(iwork)
ALLOCATE( work( lwork))
ALLOCATE(rwork(lrwork))
ALLOCATE(iwork(liwork))
call zheevd(jobz, uplo, n, matA, n, w, work, lwork, rwork, lrwork,&
            iwork, liwork, info)
call err_handler(info)

DEALLOCATE(work)
DEALLOCATE(rwork)

end subroutine zheevd_wrp
!!***

subroutine zhesv_wrp(matA, n, matB)
!! FUNCTION
!! Call to zhesv. Solve linear system A * X = B, with A a cplx*16 square
!! Hermitian matrix of size n. The matrix is destroyed during the call. B is a
!! n by n cplx*16 matrix
!!
!! INPUTS
!! matA: input matrix (cplx*16, Hermitian)
!! n: size of the system
!! matB: right hand side (cplx*16, square)
!!
!! OUTPUT
!! matA contains a factoized form of matA, see Lapack doc
!! matB contains the solutions of the system

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(INOUT), DIMENSION(n,n) :: matA, matB
!Local variables ------------------------------
CHARACTER(1) :: jobz, uplo
INTEGER :: lwork, info
COMPLEX(dp), ALLOCATABLE :: work(:)
INTEGER, ALLOCATABLE :: ipiv(:)
! *************************************************************************

! specification of the problem:
jobz = 'V' ! eigenvalues AND eigenvectors to be returned
uplo = 'U' ! consider the upper diagonal
ALLOCATE(ipiv(n))

! Query optimal work size
lwork = -1
ALLOCATE(work(1))
call zhesv(uplo,n,n,matA,n,ipiv,matB,n,work,lwork,info)
call err_handler(info)
lwork = INT(work(1))

! Actual call to zhesv
DEALLOCATE(work)
ALLOCATE(work(lwork))
call zhesv(uplo,n,n,matA,n,ipiv,matB,n,work,lwork,info)

DEALLOCATE(work)
DEALLOCATE(ipiv)

end subroutine zhesv_wrp
!!***

subroutine zgeev_vr_wrp(matA, n, w, vr)
!! FUNCTION
!! Call to zgeev. Computes the eigenvalues and right eigenvector of the cplx*16
!! square no-symmetric matrix of size n. The matrix is destroyed during the
!! call. Onnly the right eigenvectors are computed here.
!!
!! INPUTS
!! matA: input matrix (cplx*16, non-symmetric)
!! n: size of the matrix
!!
!! OUTPUT
!! matA is overwritten
!! w contains the eigenvalues
!! vr contains the right eigenvectors

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(INOUT), DIMENSION(n,n) :: matA, vr
COMPLEX(dp), INTENT(OUT), DIMENSION(n) :: w
!Local variables ------------------------------
CHARACTER(1) :: jobvl='N' ! do not compute left eigenvecs
CHARACTER(1) :: jobvr='V' ! compute right eigenvecs
INTEGER :: ldvl, lwork, lrwork, info
COMPLEX(dp), ALLOCATABLE :: work(:)
COMPLEX(dp), ALLOCATABLE :: vl(:,:)
REAL(dp), ALLOCATABLE :: rwork(:)
! *************************************************************************

! specification of the problem:
jobvl='N' ! do not compute left eigenvectors
jobvr='V' ! computes right eigenvectors
ldvl = 1 ! no space needed since we do not compute the left eigenvectors
lrwork = max(1, 2*n)
ALLOCATE(rwork(lrwork))
ALLOCATE(vl(ldvl,n))

! Query optimal work size
lwork = -1
ALLOCATE(work(1))
call zgeev(jobvl, jobvr, n, matA, n, w, vl, ldvl, vr, n, work, lwork, rwork,\
           info)
call err_handler(info)
lwork = INT(work(1))

! Actual call to zgeev
DEALLOCATE(work)
ALLOCATE(work(lwork))
call zgeev(jobvl, jobvr, n, matA, n, w, vl, ldvl, vr, n, work, lwork, rwork,\
           info)
call err_handler(info)

DEALLOCATE(work)
DEALLOCATE(rwork)
DEALLOCATE(vl)

end subroutine zgeev_vr_wrp
!!***

subroutine zggev3_vr_wrp(matA, matB, n, alpha, beta, vr)
!! FUNCTION
!! Call to zggev3. Solves a generalized eigenvalue problem A x = w B X.
!! Computes the eigenvalues and right eigenvector for the cplx*16
!! square non-symmetric matrices A and B of size n. The matrices are destroyed
!! during the call. Only the right eigenvectors are computed here.
!! Each eigenvector is scaled so the largest component has
!! abs(real part) + abs(imag. part) = 1.
!!
!! INPUTS
!! matA: input matrix (cplx*16, non-symmetric)
!! n: size of the matrix
!!
!! OUTPUT
!! matA, matB are overwritten
!! alpha divided bybeta are the eigenvalues, beta might be zero
!! vr contains the right eigenvectors

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(INOUT), DIMENSION(n,n) :: matA, matB, vr
COMPLEX(dp), INTENT(OUT), DIMENSION(n) :: alpha, beta
!Local variables ------------------------------
CHARACTER(1) :: jobvl='N' ! do not compute left eigenvecs
CHARACTER(1) :: jobvr='V' ! compute right eigenvecs
INTEGER :: ldvl, lwork, lrwork, info
COMPLEX(dp), ALLOCATABLE :: work(:)
COMPLEX(dp), ALLOCATABLE :: vl(:,:)
REAL(dp), ALLOCATABLE :: rwork(:)
! *************************************************************************

! specification of the problem:
jobvl='N' ! do not compute left eigenvectors
jobvr='V' ! computes right eigenvectors
ldvl = 1 ! no space needed since we do not compute the left eigenvectors
lrwork = max(1, 8*n)
ALLOCATE(rwork(lrwork))
ALLOCATE(vl(ldvl,n))

! Query optimal work size
lwork = -1
ALLOCATE(work(1))
call zggev3(jobvl, jobvr, n, matA, n, matB, n, alpha, beta, \
            vl, ldvl, vr, n, work, lwork, rwork, info)
call err_handler(info)
lwork = INT(work(1))

! Actual call to zggev3
DEALLOCATE(work)
ALLOCATE(work(lwork))
call zggev3(jobvl, jobvr, n, matA, n, matB, n, alpha, beta, \
            vl, ldvl, vr, n, work, lwork, rwork, info)
call err_handler(info)

DEALLOCATE(work)
DEALLOCATE(rwork)
DEALLOCATE(vl)

end subroutine zggev3_vr_wrp
!!***

subroutine zgetrf_wrp(matA, n, ipiv)
!! FUNCTION
!! Call to zgetrf. LU factorization of a cplx*16
!! square non-symmetric matrix A of size n. The LU factorization is returned
!! inside A
!!
!! INPUTS
!! matA: input matrix (cplx*16, non-symmetric)
!! n: size of the matrix
!!
!! OUTPUT
!! matA contains the LU factorization (see the Lapack doc)
!! ipiv: pivot indices (see docs)

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(INOUT), DIMENSION(n,n) :: matA
INTEGER, INTENT(OUT), DIMENSION(n) :: ipiv
!Local variables ------------------------------
INTEGER :: info
! *************************************************************************

call zgetrf(n, n, matA, n, ipiv, info)
call err_handler(info)

end subroutine zgetrf_wrp
!!***

subroutine zgetri_wrp(matA, n, ipiv)
!! FUNCTION
!! Call to zgetri. Inverses a cplx*16 square non-symmetric matrix A of size n
!! using a LU factorization as returned by zgetrf. On return A contains the
!! inverse
!!
!! INPUTS
!! matA: input matrix (cplx*16, non-symmetric)
!! n: size of the matrix
!! ipiv: pivot indices as returned by zgetrf
!!
!! OUTPUT
!! matA contains the inverse

!Arguments ------------------------------------
INTEGER, INTENT(IN) :: n
COMPLEX(dp), INTENT(INOUT), DIMENSION(n,n) :: matA
INTEGER, INTENT(IN), DIMENSION(n) :: ipiv
!Local variables ------------------------------
INTEGER :: lwork, info
COMPLEX(dp), ALLOCATABLE :: work(:)
! *************************************************************************

! Query optimal work size
lwork = -1
ALLOCATE(work(1))
call zgetri(n, matA, n, ipiv, work, lwork, info)
call err_handler(info)
lwork = INT(work(1))

! Actual call to zgetri
DEALLOCATE(work)
ALLOCATE(work(lwork))
call zgetri(n, matA, n, ipiv, work, lwork, info)
call err_handler(info)

DEALLOCATE(work)

end subroutine zgetri_wrp
!!***

END MODULE m_lapack_wrappper
