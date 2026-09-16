!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.9
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: LDT_run_amsr3_opl
! \label{LDT_amsr3_oplMod}
!
! !REVISION HISTORY:
! 16 Sep 2026: Ehsan Jalilvand; Initial Specification
!
! !INTERFACE:
subroutine LDT_run_amsr3_opl()
! !USES:
   use LDT_logMod, only: LDT_logunit
   use LDT_amsr3_oplMod, only: LDT_amsr3_oplRun

   implicit none
!
! !ARGUMENTS:
! none
!
! !DESCRIPTION:
!  This subroutine runs the runmode for resampling OPL AMSR3
!  brightness temperatures.
!EOP

   integer :: n

   n = 1
   call LDT_amsr3_oplRun(n)
   write(LDT_logunit,*) "Finished LDT AMSR3 Resampling"
end subroutine LDT_run_amsr3_opl
