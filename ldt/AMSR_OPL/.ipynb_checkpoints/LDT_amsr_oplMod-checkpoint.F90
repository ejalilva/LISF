!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: LDT_amsr_oplMod
! 
! !DESCRIPTION: 
! This module handles the run mode plugin for the 
! Operation Enhanced (9-km) AMSR soil moisture
!
! !REVISION HISTORY: 
!  14 Dec 2021: Yonghwan Kwon, Initial Specification
!  06 Feb 2023: Eric Kemp, now process subset of SMAP fields.
!  14 Feb 2023: Eric Kemp, now uses USAFSI and USAF LIS output.
!  22 Feb 2023: Eric Kemp, ensemble size now in ldt.config file.
!  01 Jul 2023: Mahdi Navari,This edit generates a separate SMAP_filelist
!                     for each LDT job based on user input.
!                     Now we can run several LDT jobs in the same directory.

! =========================
! E.J:
! check with SMAPeOPL for lower and upper case when switched from L1B and SMAP and AMSR and L1R
! watch the video with Mahdi and undrstand how the AMSR inputs are received from the config file
! EJ: why there are two 6 AM and PM mu and sigma

#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"

module LDT_amsr_oplMod

  ! Defaults
  implicit none
  private

  ! Public methods
  public :: LDT_amsr_oplInit
  public :: LDT_amsr_oplRun

  ! Public type
  type, public :: amsr_opl_dec

    character*100        :: L1Rdir, L1Rresampledir, L1Rresampledir_02, SMoutdir 
    character*100        :: LISdir, LISsnowdir
    character*100        :: TAUdir, OMEGAfile, BDfile, &
                            CLAYfile, Hfile, LCfile
    character*100        :: dailystats_ref, dailystats_lis
    character*10         :: date_curr
    integer              :: L1RresampWriteOpt, L1Rtype, AMSRfilelistSuffixNumber
    integer              :: Teffscale
    integer              :: ntimes,ngrid
    real, allocatable    :: mu_6am_ref(:), mu_6pm_ref(:) !(ngrid) ! EJ: why there are two 6 AM and PM where are they being used? why it is a vector?
    real, allocatable    :: sigma_6am_ref(:), sigma_6pm_ref(:) !(ngrid)
    real, allocatable    :: mu_6am_lis(:), mu_6pm_lis(:) !(ngrid)
    real, allocatable    :: sigma_6am_lis(:), sigma_6pm_lis(:) !(ngrid)
    integer, allocatable :: grid_col(:), grid_row(:) !(ngrid)
    real, allocatable    :: ARFS_TB_10V(:,:), ARFS_TB_10H(:,:), ARFS_TB_18H(:,:), ARFS_TB_18V(:,:), ARFS_TB_23H(:,:), ARFS_TB_23V(:,:), ARFS_TB_36H(:,:), ARFS_TB_36V(:,:), ARFS_TB_89H(:,:), ARFS_TB_89V(:,:), ARFS_LAND_WATER_FRAC(:,:) ! E.J just TB_10H but other bands should be added as well
    real                 :: SD_thold
    integer              :: num_ens ! Number of ensemble members in LIS USAF file.
    integer              :: num_tiles ! Total number of tiles in LIS USAF file.
    integer              :: ntiles_pergrid ! Number of tiles per grid point
  end type amsr_opl_dec

  type(amsr_opl_dec), public :: AMSReOPL  

contains

  subroutine LDT_amsr_oplInit()
  ! Reads Operational Enhanced AMSR-specific entries from ldt.config

    ! Imports
    use ESMF
    use LDT_coreMod, only: LDT_config
    use LDT_logMod, only: LDT_logunit, LDT_endrun, LDT_verify

    ! Defaults
    implicit none

    ! Local variables
    character(len=255) :: cfg_entry
    integer :: rc

    ! *** Get former environment variables ***

    ! Get L1Rdir
    cfg_entry = "AMSR_OPL valid date (YYYYMMDDHH):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%date_curr, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL soil moisture output directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%SMoutdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL L1R data directory:" 
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%L1Rdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL L1R data type:"    !1: NRT; 2: Historical E.J: currently there is no difference
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%L1Rtype, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL write L1R resampled output:"    !0: off; 1: on
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%L1RresampWriteOpt, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL filelist suffix number:" ! E.J: what is this in SMAP_E_OPL? do we really need it?
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%AMSRfilelistSuffixNumber, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    if(AMSReOPL%L1RresampWriteOpt.eq.1) then
       cfg_entry = "AMSR_OPL L1R resampled output directory:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%L1Rresampledir, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    endif

    ! ========= I think we can get rid of all of this with 37GHz band
    cfg_entry = "AMSR_OPL LIS soil temperature directory:" 
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%LISdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL apply soil temperature bias correction:"  !0: off; 1: on
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%Teffscale, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    if(AMSReOPL%Teffscale.eq.1) then
       cfg_entry = "AMSR_OPL reference Teff daily statistics file:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%dailystats_ref, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")

       cfg_entry = "AMSR_OPL LIS Teff daily statistics file:"
       call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
       call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%dailystats_lis, rc=rc)
       call LDT_verify(rc, trim(cfg_entry)//" not specified")
    endif
    ! ========= I think we can get rid of all of this with 37GHz band

    cfg_entry = "AMSR_OPL LIS snow directory:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%LISsnowdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL LIS ensemble size:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%num_ens, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    if (AMSReOPL%num_ens < 1) then
       write(LDT_logunit,*)'[ERR] LIS ensemble size must be at least 1!'
       write(LDT_logunit,*)'[ERR] Read in ', AMSReOPL%num_ens
       call LDT_endrun()
    end if

    cfg_entry = "AMSR_OPL LIS total number of tiles (including ensembles):"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%num_tiles, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    if (AMSReOPL%num_tiles < 1) then
       write(LDT_logunit,*) &
            '[ERR] LIS total number of tiles (including ensembles) must be'  &
            //'at least 1!'
       write(LDT_logunit,*)'[ERR] Read in ', AMSReOPL%num_tiles
       call LDT_endrun()
    end if

    cfg_entry = "AMSR_OPL LIS number of tiles per grid point:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%ntiles_pergrid, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    if (AMSReOPL%num_tiles < 1) then
       write(LDT_logunit,*) &
            '[ERR] LIS number of tiles per grid point must be at least 1!'
       write(LDT_logunit,*)'[ERR] Read in ', AMSReOPL%ntiles_pergrid
       call LDT_endrun()
    end if

    cfg_entry = "AMSR_OPL snow depth threshold:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%SD_thold, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL TAU directory:" 
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%TAUdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL OMEGA file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%OMEGAfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL soil bulk density file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%BDfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL soil clay fraction file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%CLAYfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL roughness file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%Hfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

    cfg_entry = "AMSR_OPL landcover file:"
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")
    call ESMF_ConfigGetAttribute(LDT_config, AMSReOPL%LCfile, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//" not specified")

  end subroutine LDT_amsr_oplInit


  subroutine LDT_amsr_oplRun(n)
! This calls the actual AMSR_OPL driver

! !USES:
    use esmf
    use LDT_coreMod
    use LDT_logMod

    implicit none
! !ARGUMENTS:
    integer, intent(in) :: n
!EOP

    integer, external       :: LDT_create_subdirs
    integer                 :: i, fi
    integer                 :: ftn, ierr
    character*100           :: fname
    character*100           :: amsr_l1r_filename(10) ! Why 10?
    character*8             :: yyyymmdd, yyyymmdd_01, yyyymmdd_02, yyyymmdd_03
    character*6             :: hhmmss(10)
    character*4             :: yyyy, yyyy_01, yyyy_02, yyyy_03
    character*2             :: hh, mm, dd
    character*2             :: hh_01, mm_01, dd_01
    character*2             :: hh_02, mm_02, dd_02
    character*2             :: hh_03, mm_03, dd_03
    character*2             :: tmp
    character*1             :: Orbit
    integer                 :: yr, mo, da, hr
    integer                 :: yr_pre, mo_pre, da_pre, hh_pre
    integer                 :: yr_02, mo_02, da_02, hr_02
    integer                 :: yr_03, mo_03, da_03, hr_03
    logical                 :: dir_exists, read_L1Rdata
    real                    :: teff_01(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: teff_02(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: teff_03(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: SnowDepth(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real*8                  :: TIMEsec(LDT_rc%lnc(n),LDT_rc%lnr(n))
    real                    :: UTChr(LDT_rc%lnc(n),LDT_rc%lnr(n))
    integer                 :: L1R_dir_len ! check with SMAPeopl
    integer                 :: doy_pre, doy_curr
    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: firsttime, secondtime, thirdtime, curtime, prevdaytime
    type(ESMF_TimeInterval) :: deltatime
    integer                 :: deltahr
    integer                 :: rc
    integer                 :: col, row
    external                :: readUSAFSI_amsr
    external                :: readLIS_Teff_usaf_amsr

    allocate(LDT_rc%nensem(LDT_rc%nnest))

    ! Resample AMSR L1R to L1C
    call search_AMSRL1R_files(AMSReOPL%L1Rdir,AMSReOPL%date_curr,&
                              AMSReOPL%L1Rtype, AMSReOPL%AMSRfilelistSuffixNumber)

    yyyymmdd = AMSReOPL%date_curr(1:8)
    yyyy     = AMSReOPL%date_curr(1:4)
    mm       = AMSReOPL%date_curr(5:6)
    dd       = AMSReOPL%date_curr(7:8)
    hh       = AMSReOPL%date_curr(9:10)

    if(AMSReOPL%L1RresampWriteOpt.eq.1) then
       AMSReOPL%L1Rresampledir_02 = trim(AMSReOPL%L1Rresampledir)//'/'//&
                                    trim(yyyymmdd)//'/'//trim(hh)

       ierr = LDT_create_subdirs(len_trim(AMSReOPL%L1Rresampledir_02), &
          trim(AMSReOPL%L1Rresampledir_02))
    endif

    write (tmp,'(I2.2)') AMSReOPL%AMSRfilelistSuffixNumber ! E.J: This is writing a 2 digit integer e.g., 01, 08 to tmp, but why we need this suffix for files?

    ! Reading AMSR file names one at a time for the same date and store it in amsr_L1R_filename
    ftn = LDT_getNextUnitNumber()
    open(unit=ftn, file='./AMSR_L1R_filelist_'//tmp//'.dat',&
         status='old', iostat=ierr)
    fi = 0
    do while (ierr .eq. 0)
       read (ftn, '(a)', iostat=ierr) fname
       if (ierr .ne. 0) then
          exit
       endif
       fi = fi + 1
       amsr_L1R_filename(fi) = fname
    enddo
    call LDT_releaseUnitNumber(ftn)

    L1R_dir_len = len_trim(AMSReOPL%L1Rdir)
    read_L1Rdata = .false.
    if(fi.ge.1) then
       do i=1,fi
          hhmmss(i) = trim(amsr_L1R_filename(i)(L1R_dir_len+35:L1R_dir_len+40))
          hhmmss(i+1) = trim(amsr_L1R_filename(i+1)(L1R_dir_len+35:L1R_dir_len+40))

          ! use latest version (i.e., highest version number of N**** and/or 00*) E.J: how this is reading the latest version?
          if(i == fi) then
             write (LDT_logunit,*) '[INFO] Resampling ', trim(amsr_L1R_filename(i))
             allocate(AMSReOPL%ARFS_TB_10H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_10V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_18H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_18V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_23H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_23V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_36H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_36V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_89H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_89V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_LAND_WATER_FRAC(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             ! EMK...Process subset of fields.
             call AMSR_L1R_RESAMPLE(amsr_L1R_filename(i), & 
                  AMSReOPL%L1Rdir, Orbit, TIMEsec, rc)

             if (rc == 0) then
                write (LDT_logunit,*) '[INFO] Finished resampling ', trim(amsr_L1R_filename(i))
                read_L1Rdata = .true.
             else
                deallocate(AMSReOPL%ARFS_TB_10H)
                deallocate(AMSReOPL%ARFS_TB_10V)
                deallocate(AMSReOPL%ARFS_TB_18H)
                deallocate(AMSReOPL%ARFS_TB_18V)
                deallocate(AMSReOPL%ARFS_TB_23H)
                deallocate(AMSReOPL%ARFS_TB_23V)
                deallocate(AMSReOPL%ARFS_TB_36H)
                deallocate(AMSReOPL%ARFS_TB_36V)
                deallocate(AMSReOPL%ARFS_TB_89H)
                deallocate(AMSReOPL%ARFS_TB_89V)
                deallocate(AMSReOPL%ARFS_LAND_WATER_FRAC)                
             end if
          elseif(hhmmss(i) /= hhmmss(i+1)) then
             write (LDT_logunit,*) '[INFO] Resampling ', trim(amsr_L1R_filename(i))
             allocate(AMSReOPL%ARFS_TB_10H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_10V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_18H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_18V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_23H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_23V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_36H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_36V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_89H(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_TB_89V(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             allocate(AMSReOPL%ARFS_LAND_WATER_FRAC(LDT_rc%lnc(n),LDT_rc%lnr(n)))
             !EMK Process subset of fields.
             call AMSR_L1R_RESAMPLE(amsr_L1R_filename(i), &
                  AMSReOPL%L1Rdir, Orbit, TIMEsec, rc)

             if (rc == 0) then
                write (LDT_logunit,*) '[INFO] Finished resampling ', trim(amsr_L1R_filename(i))
                read_L1Rdata = .true.
             else
                deallocate(AMSReOPL%ARFS_TB_10H)
                deallocate(AMSReOPL%ARFS_TB_10V)
                deallocate(AMSReOPL%ARFS_TB_18H)
                deallocate(AMSReOPL%ARFS_TB_18V)
                deallocate(AMSReOPL%ARFS_TB_23H)
                deallocate(AMSReOPL%ARFS_TB_23V)
                deallocate(AMSReOPL%ARFS_TB_36H)
                deallocate(AMSReOPL%ARFS_TB_36V)
                deallocate(AMSReOPL%ARFS_TB_89H)
                deallocate(AMSReOPL%ARFS_TB_89V)
                deallocate(AMSReOPL%ARFS_LAND_WATER_FRAC) 
             end if
          endif

          if(read_L1Rdata) then
  ! Get effective soil temperature (Teff) from LIS outputs

             ! use LIS outputs from previous day
             read(yyyy,*,iostat=ierr)  yr
             read(mm,*,iostat=ierr)    mo
             read(dd,*,iostat=ierr)    da
             read(hh,*,iostat=ierr)    hr

             calendar = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN, & ! define the calendar
                  name="Gregorian", &
                  rc=rc)

             ! Set current time
             call ESMF_TimeSet(curtime, yy=yr, mm=mo, dd=da, h=hr, m=0, s=0, & ! setting the current time in gregorian calendar
                  calendar=calendar, rc=rc)
             call LDT_verify(rc, '[ERR] in ESMF_TimeSet in LDT_amsr_oplRun')
             ! Go back 24 hours
             call ESMF_TimeIntervalSet(deltatime, d=1, rc=rc) ! create a time inteval = 1 day or 24 hrs
             call LDT_verify(rc, &
                  '[ERR] in ESMF_TimeIntervalSet in LDT_amsr_oplRun')
             prevdaytime = curtime - deltatime

             ! Now, find the nearest 3-hrly time (00Z, 03Z, ..., 21Z) prior
             ! to prevdaytime
             if (mod(hr, 3) == 0) then
                deltahr = 0
             else
                deltahr = hr - ((floor(real(hr)/3.))*3)
             end if
             call ESMF_TimeIntervalSet(deltatime, h=deltahr, rc=rc) ! define a new time interval based on the distance to closest 3-hourly time slot
             call LDT_verify(rc, &
                  '[ERR] in ESMF_TimeIntervalSet in LDT_amsr_oplRun')
             firsttime = prevdaytime - deltatime

             ! Now, find the next 3-hrly time (00Z, 03Z, ..., 21Z) after
             ! firsttime
             call ESMF_TimeIntervalSet(deltatime, h=3, rc=rc)
             call LDT_verify(rc, &
                  '[ERR] in ESMF_TimeIntervalSet in LDT_amsr_oplRun')
             secondtime = firsttime + deltatime

             ! Now, find the next 3-hrly time (00Z, 03Z, ..., 21Z) after
             ! secondtime
             call ESMF_TimeIntervalSet(deltatime, h=3, rc=rc)
             call LDT_verify(rc, &
                  '[ERR] in ESMF_TimeIntervalSet in LDT_amsr_oplRun')
             thirdtime = secondtime + deltatime

             ! Now, read the first time.
             call ESMF_TimeGet(firsttime, yy=yr_pre, mm=mo_pre, dd=da_pre, &
                  h=hh_pre)

             write(unit=yyyy_01, fmt='(i4.4)') yr_pre
             write(unit=mm_01, fmt='(i2.2)') mo_pre
             write(unit=dd_01, fmt='(i2.2)') da_pre
             yyyymmdd_01 = trim(yyyy_01)//trim(mm_01)//trim(dd_01)
             write(unit=hh_01, fmt='(i2.2)') hh_pre
             
             call readLIS_Teff_usaf_amsr(n, yyyymmdd_01, hh_01, Orbit, teff_01, rc)
             
             if (rc .ne. 0) then
                write(LDT_logunit,*)'[WARN] No Teff data available...'
             endif

             ! Now, read the second time.
             call ESMF_TimeGet(secondtime, yy=yr_02, mm=mo_02, dd=da_02, &
                  h=hr_02)

             write(unit=yyyy_02, fmt='(i4.4)') yr_02
             write(unit=mm_02, fmt='(i2.2)') mo_02
             write(unit=dd_02, fmt='(i2.2)') da_02
             write(unit=hh_02, fmt='(i2.2)') hr_02
             yyyymmdd_02 = trim(yyyy_02)//trim(mm_02)//trim(dd_02)
             
             call readLIS_Teff_usaf_amsr(n, yyyymmdd_02, hh_02, Orbit, teff_02, rc)
             
             if (rc .ne. 0) then
                write(LDT_logunit,*)'[WARN] No Teff data available...'
             endif

             ! Now read the third time.
             call ESMF_TimeGet(thirdtime, yy=yr_03, mm=mo_03, dd=da_03, &
                  h=hr_03)

             write(unit=yyyy_03, fmt='(i4.4)') yr_03
             write(unit=mm_03, fmt='(i2.2)') mo_03
             write(unit=dd_03, fmt='(i2.2)') da_03
             write(unit=hh_03, fmt='(i2.2)') hr_03
             yyyymmdd_03 = trim(yyyy_03)//trim(mm_03)//trim(dd_03)
             
             call readLIS_Teff_usaf_amsr(n, yyyymmdd_03, hh_03, Orbit, teff_03, rc)
             
             if (rc .ne. 0) then
                write(LDT_logunit,*)'[WARN] No Teff data available...'
             endif

             ! Now we have 3 teff for closest 3-hourly time interval in the day before AMSR retrieval (teff_01,teff_02,teff_03)
             ! Next step: Scale LIS teff to GEOS teff climatology
             ! get DOY
             call get_doy(mo_pre,da_pre,doy_pre)
             if(AMSReOPL%Teffscale.eq.1) then
                ! get getattributes
                call getattributes(AMSReOPL%dailystats_ref,&
                                   AMSReOPL%ntimes,AMSReOPL%ngrid)
                
                ! read 6-yr daily mean and std dev
                allocate(AMSReOPL%mu_6am_ref(AMSReOPL%ngrid))
                allocate(AMSReOPL%mu_6pm_ref(AMSReOPL%ngrid))
                allocate(AMSReOPL%sigma_6am_ref(AMSReOPL%ngrid))
                allocate(AMSReOPL%sigma_6pm_ref(AMSReOPL%ngrid))
                allocate(AMSReOPL%mu_6am_lis(AMSReOPL%ngrid))
                allocate(AMSReOPL%mu_6pm_lis(AMSReOPL%ngrid))
                allocate(AMSReOPL%sigma_6am_lis(AMSReOPL%ngrid))
                allocate(AMSReOPL%sigma_6pm_lis(AMSReOPL%ngrid))
                allocate(AMSReOPL%grid_col(AMSReOPL%ngrid))
                allocate(AMSReOPL%grid_row(AMSReOPL%ngrid))

                call read_DailyTeffStats_amsr(doy_pre)
                ! scale
                write (LDT_logunit,*) '[INFO] Scaling LIS effective soil temperature'
                call scale_teff_amsr(n, Orbit, teff_01, teff_02, teff_03)
                write (LDT_logunit,*) '[INFO] Finished scaling LIS effective soil temperature'

                deallocate(AMSReOPL%mu_6am_ref)
                deallocate(AMSReOPL%mu_6pm_ref)
                deallocate(AMSReOPL%sigma_6am_ref)
                deallocate(AMSReOPL%sigma_6pm_ref)
                deallocate(AMSReOPL%mu_6am_lis)
                deallocate(AMSReOPL%mu_6pm_lis)
                deallocate(AMSReOPL%sigma_6am_lis)
                deallocate(AMSReOPL%sigma_6pm_lis)
                deallocate(AMSReOPL%grid_col)
                deallocate(AMSReOPL%grid_row)
             endif
             read_L1Rdata = .false.

  ! Get snow information from LIS outputs
             call readUSAFSI_amsr(n, yyyymmdd, hh, SnowDepth, rc)
             if (rc .ne. 0) then
                write(LDT_logunit,*)'[WARN] No USAFSI data available!'
             endif

  ! Retrieve AMSR soil moisture
             ! get DOY
             call get_doy(mo,da,doy_curr)

             ! get UTC
             call get_UTC(n,TIMEsec,UTChr)

             !write(LDT_logunit,*)'EMK: UTChr = ', UTChr

             ! retrieve
             ierr = LDT_create_subdirs(len_trim(AMSReOPL%SMoutdir), &
                trim(AMSReOPL%SMoutdir))
             call ARFSSMRETRIEVAL_AMSR(amsr_L1R_filename(i), &
                  teff_01, teff_02, teff_03, &
                  SnowDepth, doy_curr, UTChr, firsttime, secondtime, thirdtime)
             deallocate(AMSReOPL%ARFS_TB_10H)
          endif
       enddo
    endif

  end subroutine LDT_amsr_oplRun

  subroutine search_AMSRL1R_files(ndir,date_curr,L1Rtype,suffix) ! EJ: in the current form there is no difference between NRT and standard product and Ascending and Descending are not separated because this seems to be a full orbit
  ! for now we are using half orbit format: GW1AM2_201207022318_175D_L1SGRTBR_2220220.h5
  ! full orbit sample from NOAA NESDIS: 20250120203202_GW1AM2_202501201817_126A_L1DLRTBR_2210210.h5


    implicit none
! !ARGUMENTS:
    character (len=*) :: ndir
    character (len=*) :: date_curr
    integer           :: L1Rtype,suffix

! !Local variables
    character*10      :: yyyymmddhh
    !character*1       :: orbit ! A or D
    character*1       :: process_kind
    character*2       :: tmp
    character*200     :: list_files

    yyyymmddhh = date_curr(1:10)

    write (tmp,'(I2.2)') suffix
    if(L1Rtype.eq.1) then   !NRT
       list_files = 'ls '//trim(ndir)//'/GW1AM2_'//&
                    trim(yyyymmddhh)// &
                    !'_*'//trim(orbit)// &
                    '*_L1SN*'//'*.h5 > AMSR_L1R_filelist_'//trim(tmp)//'.dat'
    elseif(L1Rtype.eq.2) then   !Historical
       list_files = 'ls '//trim(ndir)//'/GW1AM2_'//&
                    trim(yyyymmddhh)// &
                    !'_*'//trim(orbit)// &
                    '*_L1SG*'//'*.h5 > AMSR_L1R_filelist_'//trim(tmp)//'.dat'
    endif

    call system(trim(list_files))

  end subroutine search_AMSRL1R_files


end module LDT_amsr_oplMod
