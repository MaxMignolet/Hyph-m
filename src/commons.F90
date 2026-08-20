! small module with a few useful definitions

! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

MODULE m_commons

implicit none

PUBLIC


integer, parameter :: dp=kind(1.0d0)
real(dp), parameter :: zero=0._dp
real(dp), parameter :: one=1._dp
real(dp), parameter :: two=2._dp
real(dp), parameter :: half=0.5_dp

integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
complex(dpc), parameter :: czero = (0._dp,0._dp)
complex(dpc), parameter :: cone  = (1._dp,0._dp)
complex(dpc), parameter :: j_dpc = (0._dp,1.0_dp)

real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
real(dp), parameter :: two_pi=two*pi

real(dp), parameter :: tol5= 1.0d-05
real(dp), parameter :: tol8= 1.0d-08
real(dp), parameter :: tol12=1.0d-12
real(dp), parameter :: tol16=1.0d-16

real(dp), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit in electronic mass
real(dp), parameter :: Ha_2_cm=219474.63_dp
real(dp), parameter :: Ha_2_meV=27.211386245981d+3
real(dp), parameter :: Ha_eV=27.211386245981_dp ! 1 Hartree in eV
real(dp), parameter :: Ha_THz=6579.683920722_dp ! 1 Hartree in THz
real(dp), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
real(dp), parameter :: Bohr_Angstrom=0.529177210544_dp

character(len=1), parameter :: ch10 = char(10)
integer, parameter :: len100=100


PUBLIC :: occ_bose_einstein


contains


!!****f* m_commons/occ_bose_einstein
!! NAME
!! occ_bose_einstein
!!
!! FUNCTION
!! Computes the Bose-Einstein occupation factor for phonons
!!
!! INPUTS
!! freq: frequency of the mode (in Hartree)
!! temp: temperature
!!
!! OUTPUT
!! occ: the occupation factor
!!
!! SOURCE

function occ_bose_einstein(freq,temp) result(occ)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: freq,temp
 real(dp) :: occ

!Local variables-------------------------------

! *************************************************************************

if (temp < tol16) then
  occ = zero
else
  occ = 1/(exp(freq/(kb_HaK*temp)) - 1)
end if

end function occ_bose_einstein


END MODULE
