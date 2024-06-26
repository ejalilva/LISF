!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine setup_ALMIPIIgfrac(n)

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_vegDataMod,   only : LIS_gfrac
  use LIS_timeMgrMod

  integer, intent(in)  :: n 
  integer              :: rc

  LIS_gfrac(n)%gfracInterval  = 864000
  LIS_gfrac(n)%gfracIntervalType = "10-day" !10-day

  call LIS_registerAlarm("LIS gfrac read alarm",LIS_rc%ts, &
       LIS_gfrac(n)%gfracInterval,&
       intervalType=LIS_gfrac(n)%gfracIntervalType) 

  call ESMF_ConfigFindLabel(LIS_config,&
       "ALMIPII greenness data directory:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config, &
       LIS_gfrac(n)%gfracfile,rc=rc)
  call LIS_verify(rc,'ALMIPII greenness data directory: not defined')
end subroutine setup_ALMIPIIgfrac



