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
! !ROUTINE: scale_teff
! \label{scale_teff}
!
! !REVISION HISTORY:
!  12 JAN 2022: Yonghwan Kwon, Initial Specification
!  21 Feb 2023: Eric Kemp, added third time level
!  7 Mar 2025: Ehsan Jalilvand modified for AMSR retrieval

! !INTERFACE:
subroutine scale_teff_amsr(n, Orbit, teff_01, teff_02, teff_03)
! !USES:
  use LDT_coreMod
  use LDT_amsr_oplMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)   :: n
  character*1, intent(in) :: Orbit
  real                  :: teff_01(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                  :: teff_02(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                  :: teff_03(LDT_rc%lnc(n),LDT_rc%lnr(n))

!EOP
  integer               :: t
  integer               :: col, row
  real                  :: mu_ref, mu_lis
  real                  :: sigma_ref, sigma_lis
  real                  :: teff_01_tmp, teff_02_tmp, teff_03_tmp

  do t = 1, AMSReOPL%ngrid
     col = AMSReOPL%grid_col(t)
     row = AMSReOPL%grid_row(t)

     if (Orbit.eq."D") then
        mu_ref = AMSReOPL%mu_6am_ref(t)
        mu_lis = AMSReOPL%mu_6am_lis(t)
        sigma_ref = AMSReOPL%sigma_6am_ref(t)
        sigma_lis = AMSReOPL%sigma_6am_lis(t)
     elseif (Orbit.eq."A") then
        mu_ref = AMSReOPL%mu_6pm_ref(t)
        mu_lis = AMSReOPL%mu_6pm_lis(t)
        sigma_ref = AMSReOPL%sigma_6pm_ref(t)
        sigma_lis = AMSReOPL%sigma_6pm_lis(t)
     endif

     if (teff_01(col,row).ne.-9999.0) then
        if (mu_ref.gt.0.and.&
             mu_lis.gt.0.and.&
             sigma_ref.gt.0.and.&
             sigma_lis.gt.0) then

           teff_01_tmp = (teff_01(col,row) - mu_lis) * &
                         sigma_ref/sigma_lis + mu_ref

           teff_01(col,row) = teff_01_tmp
        endif
     endif

     if (teff_02(col,row).ne.-9999.0) then
        if (mu_ref.gt.0.and.&
             mu_lis.gt.0.and.&
             sigma_ref.gt.0.and.&
             sigma_lis.gt.0) then

             teff_02_tmp = (teff_02(col,row) - mu_lis) * &
                         sigma_ref/sigma_lis + mu_ref

             teff_02(col,row) = teff_02_tmp
        endif
     endif

     if (teff_03(col,row).ne.-9999.0) then
        if (mu_ref.gt.0.and.&
             mu_lis.gt.0.and.&
             sigma_ref.gt.0.and.&
             sigma_lis.gt.0) then

             teff_03_tmp = (teff_03(col,row) - mu_lis) * &
                         sigma_ref/sigma_lis + mu_ref

             teff_03(col,row) = teff_03_tmp
        endif
     endif

  enddo

end subroutine scale_teff_amsr
