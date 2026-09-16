!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: LDT_amsr3_oplMod
!
! !DESCRIPTION:
!  Run-mode plugin that resamples AMSR3 L1R brightness temperatures onto
!  the Air Force (ARFS) grid and writes them as netCDF. No retrieval.
!
!  ldt.config entries (AMSR3_OPL ...):
!    required : valid date (YYYYMMDDHH), L1R data directory,
!               resampled output directory
!    optional : write resampled output (0/1, default 1)
!               filter snow and precip footprints (0/1, default 1)
!               search radius km by band (3 values, default 20 15 10)
!
!  The per-job file list is AMSR3_L1R_filelist_<YYYYMMDDHH>.dat, keyed on
!  the valid date exactly as WSF_OPL does, so concurrent LDT jobs in one
!  directory never overwrite each other's list (LISF PR #1799).
!
! !REVISION HISTORY:
! 16 Sep 2026: Ehsan Jalilvand; Initial Specification
!-------------------------------------------------------------------------

#include "LDT_misc.h"

module LDT_amsr3_oplMod

  use LDT_constantsMod, only: LDT_CONST_PATH_LEN

  implicit none
  private

  public :: LDT_amsr3_oplInit
  public :: LDT_amsr3_oplRun

  type, public :: amsr3_opl_dec
     character(len=LDT_CONST_PATH_LEN) :: L1Rdir, outdir, outdir_date
     character*10 :: date_curr
     integer      :: write_opt
     logical      :: filter_snow_precip
     real         :: radius_by_band(3)     ! km for f<=11 GHz, 11<f<=30 GHz, f>30 GHz
  end type amsr3_opl_dec

  type(amsr3_opl_dec), public :: AMSR3eOPL

contains

  subroutine LDT_amsr3_oplInit()
    use ESMF
    use LDT_coreMod, only: LDT_config
    use LDT_logMod,  only: LDT_logunit, LDT_verify
    implicit none
    character(len=255) :: cfg_entry
    integer :: rc, ival

    write(LDT_logunit,*) '[INFO] Initializing AMSR3 L1R resampling'

    cfg_entry = 'AMSR3_OPL valid date (YYYYMMDDHH):'
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//' not specified')
    call ESMF_ConfigGetAttribute(LDT_config, AMSR3eOPL%date_curr, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//' not specified')

    cfg_entry = 'AMSR3_OPL L1R data directory:'
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//' not specified')
    call ESMF_ConfigGetAttribute(LDT_config, AMSR3eOPL%L1Rdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//' not specified')

    cfg_entry = 'AMSR3_OPL resampled output directory:'
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//' not specified')
    call ESMF_ConfigGetAttribute(LDT_config, AMSR3eOPL%outdir, rc=rc)
    call LDT_verify(rc, trim(cfg_entry)//' not specified')

    ! ---- optional entries with defaults (WSF_OPL pattern) ----
    AMSR3eOPL%write_opt = 1
    cfg_entry = 'AMSR3_OPL write resampled output:'
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    if (rc == 0) call ESMF_ConfigGetAttribute(LDT_config, AMSR3eOPL%write_opt, rc=rc)

    ival = 1
    cfg_entry = 'AMSR3_OPL filter snow and precip footprints:'
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    if (rc == 0) call ESMF_ConfigGetAttribute(LDT_config, ival, rc=rc)
    AMSR3eOPL%filter_snow_precip = (ival == 1)

    AMSR3eOPL%radius_by_band = (/ 20.0, 15.0, 10.0 /)
    cfg_entry = 'AMSR3_OPL search radius km by band:'
    call ESMF_ConfigFindLabel(LDT_config, trim(cfg_entry), rc=rc)
    if (rc == 0) call ESMF_ConfigGetAttribute(LDT_config, AMSR3eOPL%radius_by_band, count=3, rc=rc)

    write(LDT_logunit,*) '[INFO] AMSR3_OPL date=', trim(AMSR3eOPL%date_curr), &
         ' filter_snow_precip=', AMSR3eOPL%filter_snow_precip, &
         ' radius_by_band=', AMSR3eOPL%radius_by_band
  end subroutine LDT_amsr3_oplInit


  subroutine LDT_amsr3_oplRun(n)
    use LDT_coreMod
    use LDT_logMod
    implicit none
    integer, intent(in) :: n
    integer, external :: LDT_create_subdirs
    integer :: ftn, ierr, rc, nfiles
    character*10 :: tmp
    character(len=LDT_CONST_PATH_LEN) :: fname

    write(LDT_logunit,*) '[INFO] Starting AMSR3 L1R resampling for ', AMSR3eOPL%date_curr

    ! file list keyed on the valid date: no collisions between concurrent jobs
    tmp = trim(AMSR3eOPL%date_curr)
    call search_AMSR3L1R_files(AMSR3eOPL%L1Rdir, AMSR3eOPL%date_curr, tmp)

    AMSR3eOPL%outdir_date = trim(AMSR3eOPL%outdir)//'/'//AMSR3eOPL%date_curr(1:8)
    if (AMSR3eOPL%write_opt == 1) &
         ierr = LDT_create_subdirs(len_trim(AMSR3eOPL%outdir_date), trim(AMSR3eOPL%outdir_date))

    ftn = LDT_getNextUnitNumber()
    open(ftn, file='AMSR3_L1R_filelist_'//trim(tmp)//'.dat', status='old', iostat=ierr)
    if (ierr /= 0) then
       write(LDT_logunit,*) '[ERR] Cannot open AMSR3_L1R_filelist_'//trim(tmp)//'.dat'
       return
    endif

    ! every granule of the day; the driver splits each into asc / desc passes
    nfiles = 0
    do while (ierr == 0)
       read(ftn, '(a)', iostat=ierr) fname
       if (ierr /= 0) exit
       if (len_trim(fname) == 0) cycle
       if (index(fname, 'No such file') > 0 .or. index(fname, 'cannot access') > 0) cycle   ! ls errors (2>&1)
       nfiles = nfiles + 1
       write(LDT_logunit,*) '[INFO] Resampling ', trim(fname)
       call AMSR3_L1R_RESAMPLE(trim(fname), trim(AMSR3eOPL%outdir_date), AMSR3eOPL%write_opt, &
                               AMSR3eOPL%filter_snow_precip, AMSR3eOPL%radius_by_band, rc)
       if (rc /= 0) write(LDT_logunit,*) '[WARN] resample failed for ', trim(fname)
    end do
    call LDT_releaseUnitNumber(ftn)
    if (nfiles == 0) write(LDT_logunit,*) '[WARN] no GGWAM3 L1R files found for ', AMSR3eOPL%date_curr(1:8)
    write(LDT_logunit,*) '[INFO] AMSR3 L1R resampling done: ', nfiles, ' granules'
  end subroutine LDT_amsr3_oplRun


  subroutine search_AMSR3L1R_files(ndir, date_curr, suffix)
    use LDT_logMod, only: LDT_logunit
    implicit none
    character(len=*) :: ndir, date_curr, suffix
    character(len=LDT_CONST_PATH_LEN) :: list_files, search_pattern
    external :: system

    ! standard L1R granule: GGWAM3_YYYYMMDDHHMM<A|D>NNN_L1RTBR....nc ; match by day
    search_pattern = trim(ndir)//'/GGWAM3_'//date_curr(1:8)//'*_L1RTBR*.nc'
    list_files = 'ls '//trim(search_pattern)//' > AMSR3_L1R_filelist_'//trim(suffix)//'.dat 2>&1'

    write(LDT_logunit,*) '[INFO] Searching AMSR3 L1R files: ', trim(search_pattern)
    call system(trim(list_files))
  end subroutine search_AMSR3L1R_files

end module LDT_amsr3_oplMod
