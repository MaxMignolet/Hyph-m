! Copyright 2026 M. Mignolet
! SPDX-License-Identifier: AGPL-3.0-or-later

MODULE m_readInput

implicit none


private

public :: readInput

CONTAINS


subroutine readInput(crystName)
!! FUNCTION
!! Read the input file 'input' and returns the name of the crystal (CrI3 or
!! CrBr3)
!!
!! INPUTS
!!
!! OUTPUT
!! crystName: name of the crystal inn the input

!Arguments ------------------------------------
character(len=100),intent(out) :: crystName
!Local variables ------------------------------
character(len=100) :: fname='input'
integer :: unit
! *************************************************************************

open(newunit=unit, file=fname)
read(unit,'(a)') crystName

return

end subroutine

END MODULE
