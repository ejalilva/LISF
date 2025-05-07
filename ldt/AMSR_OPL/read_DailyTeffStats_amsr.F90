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
! !ROUTINE: get_doy_amsr
! \label{get_doy_amsr}
!
! !INTERFACE:
subroutine get_doy_amsr(mo,da,doy)
! 
! !USES:

  implicit none

! !ARGUMENTS:
  integer    :: mo, da, doy
  integer    :: imo

!EOP

doy = 0
do imo = 1,mo
   if(imo.lt.mo) then
      if(imo.eq.1.or.&
       imo.eq.3.or.&
       imo.eq.5.or.&
       imo.eq.7.or.&
       imo.eq.8.or.&
       imo.eq.10.or.&
       imo.eq.12) then
         doy = doy + 31
      elseif(imo.eq.2) then
         doy = doy + 28
      else
         doy = doy + 30
      endif
   elseif(imo.eq.mo) then
      if(imo.eq.2) then
         if(da.eq.29) then
            da = 28
         endif
      endif
      doy = doy + da
   endif
enddo

end subroutine get_doy_amsr

!BOP
! 
! !ROUTINE: get_UTC_amsr
! \label{get_UTC_amsr}
!
! !INTERFACE:
subroutine get_UTC_amsr(n,TIMEsec,UTChr)
! 
! !USES:
  use LDT_coreMod
  use LDT_logMod, only: LDT_logunit

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  real*8              :: TIMEsec(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                :: UTChr(LDT_rc%lnc(n),LDT_rc%lnr(n))

!EOP
  integer             :: ilat, ilon, imo, ida
  real*8              :: TIMEday
  real                :: TIMEhr
  real                :: UTCyr, UTCmo, UTCda
  integer             :: count_yr, dayremove

  do ilat=1,LDT_rc%lnr(n)
     do ilon=1,LDT_rc%lnc(n)

        if (TIMEsec(ilon,ilat).gt.0) then
           !write(LDT_logunit,*) 'EMK: ilon,ilat, TIMEsec= ', &
           !     ilon, ilat, TIMEsec(ilon,ilat)

           TIMEday = TIMEsec(ilon,ilat)/DBLE(60*60*24)
           !write(*,*) 'TIMEday= ', TIMEday

           count_yr = 0
           do while (TIMEday.ge.366)
              if(mod(count_yr,4).eq.0) then
                 TIMEday = TIMEday - 366
              else
                 TIMEday = TIMEday - 365
              endif
              count_yr = count_yr + 1
           enddo
           UTCyr = 1993 + count_yr ! E.J: AMSR Scanning time data is stored as TAI93 in AMSR2 product. TAI93 is the elapsed second time which includes the leap second from January 1st, 1993.

           !write(*,*) 'UTCyr= ', UTCyr
           !write(*,*) 'TIMEday= ', TIMEday          
 
           imo = 1
           do while (TIMEday.gt.0)
              if(imo.eq.1.or.&
               imo.eq.3.or.&
               imo.eq.5.or.&
               imo.eq.7.or.&
               imo.eq.8.or.&
               imo.eq.10.or.&
               imo.eq.12) then
                 TIMEday = TIMEday - 31
                 dayremove = 31
              elseif(imo.eq.2) then
                 if(mod(UTCyr,4.).eq.0) then
                    TIMEday = TIMEday - 29
                    dayremove = 29
                 else
                    TIMEday = TIMEday - 28
                    dayremove = 28
                 endif
              else
                 TIMEday = TIMEday - 30
                 dayremove = 30
              endif

              imo = imo + 1
           enddo  
          
           if(TIMEday.lt.0) then
              TIMEday = TIMEday + dayremove
              imo = imo - 1
           endif
           UTCmo = imo
           UTCda = ceiling(TIMEday)
           TIMEhr = (TIMEday - UTCda + 1)*24
           UTChr(ilon,ilat) = TIMEhr + 12  !reference time: Jan.1, 2000 at 12pm.

           if(UTChr(ilon,ilat).gt.24) then
              UTCda = UTCda + 1
              UTChr(ilon,ilat) = UTChr(ilon,ilat) - 24
           endif

           !write(LDT_logunit,*) 'EMK: ilon,ilat,TIMEsec,UTChr= ', &
           !     ilon, ilat, TIMEsec(ilon,ilat), UTChr(ilon,ilat)

        else
           UTChr(ilon,ilat) = LDT_rc%udef 
        endif

     enddo
  enddo

write(LDT_logunit,*) '[DEBUG] UTChr statistics:'
write(LDT_logunit,*) '  Min UTChr:', MINVAL(UTChr, MASK=(UTChr > -9990))
write(LDT_logunit,*) '  Max UTChr:', MAXVAL(UTChr, MASK=(UTChr > -9990))
write(LDT_logunit,*) '  Count of negative UTChr:', COUNT(UTChr < 0)
write(LDT_logunit,*) '  Count of undefined UTChr:', COUNT(UTChr < -9990)
    
end subroutine get_UTC_amsr
