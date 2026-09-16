!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: AMSR3_L1R_RESAMPLE
!
! DESCRIPTION: Resample one AMSR3 L1R granule onto the Air Force (ARFS)
!   grid and write the result as netCDF.
!   - all 14 channels are read at native FOV (3D, channel table in reader)
!   - the granule is split into ascending / descending by the nadir
!     latitude trend, and one file is written per pass
!   - per-channel radius comes from radius_by_band; snow/precip filtering
!     is a switch; both are meant to be set from ldt.config by the mod
!   - output variable names are derived from the reader's channel names:
!     Tb_FOV06Ch06V_P89o -> TB_06V, FRFI_06V, FCAU_06V, QOR_06V
!-------------------------------------------------------------------------

subroutine AMSR3_L1R_RESAMPLE(AMSR3FILE, outdir, write_opt, filter_snow_precip, radius_by_band, rc)

  USE TOOLSUBS_AMSR3
  USE invdist_l1r2arfs_AMSR3
  USE LDT_logMod
  USE ESMF
  USE netcdf
  IMPLICIT NONE

  character(*), intent(in)  :: AMSR3FILE, outdir
  integer,      intent(in)  :: write_opt              ! 1 = write netCDF
  logical,      intent(in)  :: filter_snow_precip
  real*4,       intent(in)  :: radius_by_band(3)      ! km for f<=11 GHz, 11<f<=30 GHz, f>30 GHz
  integer,      intent(out) :: rc

  ! ---- swath (from reader) ----
  real*4,    allocatable :: tb(:,:,:), lat(:,:), lon(:,:)
  integer*2, allocatable :: tbq(:,:,:), scan_qual(:)
  real*8,    allocatable :: tb_time(:)
  integer*4, allocatable :: land_frac(:,:), snow(:,:), precip(:,:)
  character(len=20) :: ch_name(NCH) ; real*4 :: ch_freq(NCH) ; character(len=1) :: ch_pol(NCH)
  integer :: nfov, nscan, nchan, ierr

  ! ---- grid (from resampler) ----
  real*8,    allocatable :: ARFS_LAT(:), ARFS_LON(:), arfs_tim(:,:)
  real*4,    allocatable :: arfs_tb(:,:,:), arfs_frfi(:,:,:), arfs_fcau(:,:,:)
  real*4,    allocatable :: arfs_fsnow(:,:), arfs_fprec(:,:), arfs_land(:,:)
  integer*2, allocatable :: arfs_qor(:,:,:)
  integer*4, allocatable :: arfs_nsamp(:,:)

  integer, parameter :: arfs_nrow_lat = 1920, arfs_mcol_lon = 2560   ! ARFS grid size
  real*4  :: radius_km(NCH)
  logical, allocatable :: asc(:), use_scan(:)
  integer :: k, ii, ipass, nuse
  character(len=1) :: pass
  character(len=200) :: basename
  integer :: ncid, dlon, dlat                     ! shared by write_pass and def2d

  rc = 0
  write(LDT_logunit,*)'[INFO] AMSR3_L1R_RESAMPLE: ', trim(AMSR3FILE)

  ! ---- ARFS grid, built inline as WSF_OPL does (no shared grid modules) ----
  ! 2560 x 1920 cells, 0.140625 x 0.09375 deg, cell centres; lat runs south -> north
  ! (LIS convention used by AMSR2). WSF writes north -> south; flip here if the
  ! downstream reader expects the WSF orientation.
  allocate(ARFS_LAT(arfs_nrow_lat), ARFS_LON(arfs_mcol_lon))
  do k = 1, arfs_nrow_lat ; ARFS_LAT(k) = -90.d0  + (k-0.5d0)*0.09375d0  ; end do
  do k = 1, arfs_mcol_lon ; ARFS_LON(k) = -180.d0 + (k-0.5d0)*0.140625d0 ; end do

  ! ---- read the granule ----
  call get_amsr3_l1r(AMSR3FILE, tb, tbq, ch_name, ch_freq, ch_pol, lat, lon, tb_time, &
                     land_frac, snow, precip, scan_qual, nfov, nscan, nchan, ierr)
  if (ierr /= 0 .or. nscan == 0) then
     write(LDT_logunit,*)'[ERR] get_amsr3_l1r failed for ', trim(AMSR3FILE) ; rc = 1 ; return
  endif

  ! ---- per-channel radius from band ----
  do k = 1, nchan
     if (ch_freq(k) <= 11.0) then ; radius_km(k) = radius_by_band(1)
     elseif (ch_freq(k) <= 30.0) then ; radius_km(k) = radius_by_band(2)
     else ; radius_km(k) = radius_by_band(3) ; endif
  end do

  ! ---- ascending / descending per scan from the nadir latitude trend ----
  allocate(asc(nscan), use_scan(nscan))
  asc(1:nscan-1) = lat(nfov/2+1, 2:nscan) > lat(nfov/2+1, 1:nscan-1)
  asc(nscan)     = asc(nscan-1)
  write(LDT_logunit,*)'[INFO] scans ascending=', count(asc), ' descending=', count(.not. asc)

  basename = AMSR3FILE(index(AMSR3FILE, '/', back=.true.)+1:)

  do ipass = 1, 2
     if (ipass == 1) then ; pass = 'A' ; use_scan = asc
     else                 ; pass = 'D' ; use_scan = .not. asc ; endif
     nuse = count(use_scan)
     if (nuse < 50) then
        write(LDT_logunit,*)'[INFO] pass ', pass, ': only ', nuse, ' scans, skipped' ; cycle
     endif

     call L1RTB2ARFS_INVDIS_AMSR3(tb, tbq, snow, precip, scan_qual, use_scan, land_frac, tb_time, &
          lat, lon, nfov, nscan, nchan, ARFS_LAT, ARFS_LON, radius_km, filter_snow_precip, &
          arfs_tb, arfs_frfi, arfs_fcau, arfs_qor, arfs_fsnow, arfs_fprec, arfs_tim, arfs_land, arfs_nsamp)

     if (write_opt == 1) call write_pass(pass)
  end do

CONTAINS

  ! Tb_FOV06Ch06V_P89o -> '06V' ; Tb_FOV10Ch10uH_P89o -> '10uH'
  function chtag(name) result(tag)
    character(*), intent(in) :: name ; character(len=8) :: tag
    integer :: p1, p2
    p1 = index(name, 'Ch') + 2 ; p2 = index(name, '_P89o') - 1
    tag = name(p1:p2)
  end function chtag

  subroutine write_pass(pass)
    character(len=1), intent(in) :: pass
    character(len=300) :: fname
    integer :: dtim, vlon, vlat, vtim
    integer :: vtb(NCH), vrfi(NCH), vcau(NCH), vqor(NCH), vland, vsnow, vprec, vns, vstim
    integer :: yr, mo, dy, hr, mi, ios
    type(ESMF_Time) :: t0, t1 ; type(ESMF_TimeInterval) :: dt ; real*8 :: tsec

    fname = trim(outdir)//'/AMSR3_L1R_resampled_'//basename(8:19)//'_'//pass//'.nc'
    call LDT_verify(nf90_create(trim(fname), NF90_NETCDF4, ncid), '[ERR] nf90_create '//trim(fname))
    call LDT_verify(nf90_def_dim(ncid, 'lon',  arfs_mcol_lon, dlon), 'def_dim lon')
    call LDT_verify(nf90_def_dim(ncid, 'lat',  arfs_nrow_lat, dlat), 'def_dim lat')
    call LDT_verify(nf90_def_dim(ncid, 'time', 1,             dtim), 'def_dim time')

    call LDT_verify(nf90_def_var(ncid, 'lon',  NF90_FLOAT,  dlon, vlon), 'def lon')
    call LDT_verify(nf90_def_var(ncid, 'lat',  NF90_FLOAT,  dlat, vlat), 'def lat')
    call LDT_verify(nf90_def_var(ncid, 'time', NF90_DOUBLE, dtim, vtim), 'def time')
    call LDT_verify(nf90_put_att(ncid, vlon, 'units', 'degrees_east'),  'att lon')
    call LDT_verify(nf90_put_att(ncid, vlat, 'units', 'degrees_north'), 'att lat')
    call LDT_verify(nf90_put_att(ncid, vtim, 'units', 'seconds since 1970-01-01 00:00:00'), 'att time')

    do k = 1, nchan                                              ! four layers per channel
       vtb(k)  = def2d('TB_'  //chtag(ch_name(k)), NF90_FLOAT, 'K',  'brightness temperature, '//trim(ch_name(k)))
       vrfi(k) = def2d('FRFI_'//chtag(ch_name(k)), NF90_FLOAT, '1',  'weighted share of footprints with RFI occurred')
       vcau(k) = def2d('FCAU_'//chtag(ch_name(k)), NF90_FLOAT, '1',  'weighted share of footprints with resampling caution')
       vqor(k) = def2d('QOR_' //chtag(ch_name(k)), NF90_SHORT, '1',  'bitwise OR of contributing L1R quality bytes')
    end do
    vland = def2d('LAND_FRAC', NF90_FLOAT, '%', 'land area percent, FOV06')
    vsnow = def2d('FSNOW',     NF90_FLOAT, '1', 'weighted share of footprints flagged snow')
    vprec = def2d('FPRECIP',   NF90_FLOAT, '1', 'weighted share of footprints flagged precipitation')
    vns   = def2d('NSAMP',     NF90_INT,   '1', 'footprints within the largest search radius')
    vstim = def2d('SCAN_TIME', NF90_DOUBLE,'seconds since 1993-01-01 00:00:00 (TAI)', 'weighted mean scan time')
    call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'source_file', trim(basename)), 'att source')
    call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'orbit_pass',  pass), 'att pass')
    call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'search_radius_km_by_band', radius_by_band), 'att radius')
    call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'filter_snow_precip', merge(1,0,filter_snow_precip)), 'att filter')
    call LDT_verify(nf90_enddef(ncid), 'enddef')

    call LDT_verify(nf90_put_var(ncid, vlon, real(ARFS_LON)), 'put lon')
    call LDT_verify(nf90_put_var(ncid, vlat, real(ARFS_LAT)), 'put lat')

    ! scalar time: granule timestamp from the filename (YYYYMMDDHHMM at 8:19), seconds since 1970
    read(basename(8:11),'(i4)',iostat=ios) yr ; read(basename(12:13),'(i2)',iostat=ios) mo
    read(basename(14:15),'(i2)',iostat=ios) dy ; read(basename(16:17),'(i2)',iostat=ios) hr
    read(basename(18:19),'(i2)',iostat=ios) mi
    tsec = -9999.d0
    if (ios == 0) then
       call ESMF_TimeSet(t0, yy=1970, mm=1, dd=1, h=0, m=0, s=0, rc=ios)
       call ESMF_TimeSet(t1, yy=yr, mm=mo, dd=dy, h=hr, m=mi, s=0, rc=ios)
       dt = t1 - t0 ; call ESMF_TimeIntervalGet(dt, s_r8=tsec, rc=ios)
    endif
    call LDT_verify(nf90_put_var(ncid, vtim, tsec), 'put time')

    do k = 1, nchan
       call LDT_verify(nf90_put_var(ncid, vtb(k),  arfs_tb(:,:,k)),   'put TB')
       call LDT_verify(nf90_put_var(ncid, vrfi(k), arfs_frfi(:,:,k)), 'put FRFI')
       call LDT_verify(nf90_put_var(ncid, vcau(k), arfs_fcau(:,:,k)), 'put FCAU')
       call LDT_verify(nf90_put_var(ncid, vqor(k), arfs_qor(:,:,k)),  'put QOR')
    end do
    call LDT_verify(nf90_put_var(ncid, vland, arfs_land),  'put LAND_FRAC')
    call LDT_verify(nf90_put_var(ncid, vsnow, arfs_fsnow), 'put FSNOW')
    call LDT_verify(nf90_put_var(ncid, vprec, arfs_fprec), 'put FPRECIP')
    call LDT_verify(nf90_put_var(ncid, vns,   arfs_nsamp), 'put NSAMP')
    call LDT_verify(nf90_put_var(ncid, vstim, arfs_tim),   'put SCAN_TIME')
    call LDT_verify(nf90_close(ncid), 'close')
    write(LDT_logunit,*)'[INFO] wrote ', trim(fname)
  end subroutine write_pass

  ! define one (lon,lat) variable with fill, units and long_name; returns its varid
  integer function def2d(name, xtype, units, long_name)
    character(*), intent(in) :: name, units, long_name ; integer, intent(in) :: xtype
    call LDT_verify(nf90_def_var(ncid, trim(name), xtype, [dlon, dlat], def2d), 'def '//trim(name))
    select case (xtype)
    case (NF90_FLOAT)  ; call LDT_verify(nf90_def_var_fill(ncid, def2d, 0, -9999.0),  'fill '//trim(name))
    case (NF90_DOUBLE) ; call LDT_verify(nf90_def_var_fill(ncid, def2d, 0, -9999.d0), 'fill '//trim(name))
    case (NF90_SHORT)  ; call LDT_verify(nf90_def_var_fill(ncid, def2d, 0, -1_2),     'fill '//trim(name))
    case (NF90_INT)    ; call LDT_verify(nf90_def_var_fill(ncid, def2d, 0, 0),        'fill '//trim(name))
    end select
    call LDT_verify(nf90_put_att(ncid, def2d, 'units', units), 'att units '//trim(name))
    call LDT_verify(nf90_put_att(ncid, def2d, 'long_name', long_name), 'att long_name '//trim(name))
  end function def2d

end subroutine AMSR3_L1R_RESAMPLE
