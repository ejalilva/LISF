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
! !ROUTINE: LDT_init_amsr3_opl
! \label{LDT_amsr3_oplMod}
!
! !REVISION HISTORY:
! 16 Sep 2026: Ehsan Jalilvand; Initial Specification
!
! !INTERFACE:
subroutine LDT_init_amsr3_opl()
! !USES:
   use LDT_domainMod, only: LDT_setDomainSpecs
   use LDT_logMod, only: LDT_logunit
   use LDT_paramProcMod, only: LDT_paramProcConfig
   use LDT_amsr3_oplMod, only: LDT_amsr3_oplInit

   implicit none
!
! !ARGUMENTS:
! none
!
! !DESCRIPTION:
!  This subroutine initializes the runmode for resampling OPL AMSR3
!  brightness temperatures.
!EOP

   write(LDT_logunit,*) "Start of AMSR3 L1R Resampling"
   call LDT_setDomainSpecs()
   call LDT_paramProcConfig()
   call LDT_amsr3_oplInit()
   flush(LDT_logunit)
end subroutine LDT_init_amsr3_opl
