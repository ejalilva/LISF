!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: scale_teff_amsr
! \label{scale_teff_amsr}
!
! !REVISION HISTORY:
!  26 Feb 2025: Ehsan Jalilvand, modified for AMSR retrieval - simplified version

! !INTERFACE:
subroutine scale_teff_amsr(n, Orbit, teff_01, teff_02, teff_03)
! !USES:
  use LDT_coreMod
  
  implicit none
! !ARGUMENTS:
  integer, intent(in)     :: n
  character*1, intent(in) :: Orbit
  real                    :: teff_01(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                    :: teff_02(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                    :: teff_03(LDT_rc%lnc(n),LDT_rc%lnr(n))

!EOP
  ! In AMSR implementation, we don't need to perform scaling
  ! This simplified subroutine is kept for compatibility but doesn't modify data
  
end subroutine scale_teff_amsr