!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: AMSR_L1R_SUBSET
!
! REVISION HISTORY:
!  22 Oct 2021: P.W.LIU; Initial SMAP implemetation
!  18 Dec 2021: Yonghwan Kwon; modified for LDT
!  09 Feb 2023: Eric Kemp; now processes subset of fields, no output of
!    data to separate binary files.
!  19 Apr 2024: Eric Kemp; changed pixel times to REAL*8.
!  26 Feb 2025: Ehsan Jalilvand; AMSR2 implementation
!
! DESCRIPTION: RESAMPLE AMSR_L1R TB TO AIR FORCE GRID
! INPUT : AMSR - L1R Brightness Temperature
! OUTPUT: AMSRTB_ARFSGRIDE_ddmmyyy.dat
! NOTES : Inverse Distance Squared with 0.4 deg serching window
!-------------------------------------------------------------------------

subroutine AMSR_L1R_RESAMPLE(AMSRFILE,L1R_dir,Orbit,ARFS_TIME,rc)

  USE VARIABLES
  USE DATADOMAIN
  USE FUNCTIONS
  USE TOOLSUBS_AMSR
  USE invdist_l1r2arfs
  USE LDT_logMod
  USE LDT_amsr_oplMod
  USE ESMF
  USE netcdf
  
  IMPLICIT NONE

  INTEGER :: i, j, nrow, mcol
  CHARACTER (len=100) :: AMSRFILE
  character (len=100) :: L1R_dir
  character (len=20)  :: variable_name(13)
  character (len=100) :: resample_filename(12)
  character (len=1)   :: Orbit ! E.J: orbit is one of the outputs extracted from the filename (D: Descending & A: Ascending)
  integer             :: var_i
  integer             :: L1R_dir_len,L1R_fname_len
  integer :: ierr
  integer :: rc

  REAL*8,DIMENSION(:),ALLOCATABLE :: TIME_L1R
  REAL*4,DIMENSION(:,:),ALLOCATABLE :: TB_10H, TB_10V, TB_18H, TB_18V, TB_23H, TB_23V, TB_36H, TB_36V, TB_89H, TB_89V 
  REAL*4,DIMENSION(:,:),ALLOCATABLE :: LAT_L1R, LON_L1R, LAT89, LON89
  INTEGER*2,DIMENSION(:,:),ALLOCATABLE :: RFI_FLAG
  INTEGER*4,DIMENSION(:,:),ALLOCATABLE :: LAND_WATER_FRAC
  INTEGER*4,DIMENSION(:,:),ALLOCATABLE :: SNOW, PRECIP
  INTEGER*1,DIMENSION(:,:),ALLOCATABLE :: QUALITY_FLAG               ! Footprint level
  INTEGER :: nrow89,ncol89 ! m89 & n89 are for the 89GHz band for which lat and lon are provided

  REAL*8,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
  INTEGER*4,DIMENSION(2560,1920) :: ARFS_SAMPLE_V, ARFS_SAMPLE_H 
  REAL*8,DIMENSION(2560,1920) :: ARFS_TIME
  REAL*4,DIMENSION(2560,1920) :: ARFS_TB_10V, ARFS_TB_10H, ARFS_TB_18H, ARFS_TB_18V,ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, ARFS_TB_89H, ARFS_TB_89V, ARFS_LAND_WATER_FRAC ! This is instead of ARFS_COR_TBV
  INTEGER*1,DIMENSION(2560,1920) :: ARFS_QUALITY_FLAG
  !INTEGER*4,DIMENSION(2560,1920),ALLOCATABLE :: ARFS_RFI_FLAG

  REAL :: T1, T2
  character(len=200) :: netcdf_filename
  character(len=100) ::  basename, extracted_part
  character(len=12) :: datetime_str  
  integer :: filename_start_pos
  integer :: year, month, day, hour, minute
  type(ESMF_Time) :: file_time, reference_time
  type(ESMF_TimeInterval) :: time_diff
  real*8 :: time_seconds
  integer :: rc_time
  integer :: ncid, time_dimid, lat_dimid, lon_dimid
  integer :: time_varid, lat_varid, lon_varid
  integer :: tb_10h_varid, tb_10v_varid, tb_18h_varid, tb_18v_varid
  integer :: tb_23h_varid, tb_23v_varid, tb_36h_varid, tb_36v_varid  
  integer :: tb_89h_varid, tb_89v_varid, lwf_varid, qf_varid
  integer :: iret
  real, allocatable :: lats(:), lons(:)

  rc = 0
  ! Extra logging for debug
  write(LDT_logunit,*) '[INFO] Starting AMSR_L1R_RESAMPLE for file:', trim(AMSRFILE)
  
  CALL ARFS_GEO
  ALLOCATE(ARFS_LAT(arfs_nrow_lat),ARFS_LON(arfs_mcol_lon)) 
  ARFS_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_space)
  ARFS_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_space)

  write(LDT_logunit,*) '[INFO] Calling get_amsr_l1r to read:', trim(AMSRFILE)

  CALL get_amsr_l1r (AMSRFILE,TIME_L1R, &
          TB_10V, TB_10H, TB_18V, TB_18H, &
          TB_23V, TB_23H, TB_36V, TB_36H, &
          TB_89V, TB_89H, &
          LAT_L1R, LON_L1R, LAT89, LON89, &
          LAND_WATER_FRAC, SNOW, PRECIP, &
          RFI_FLAG, QUALITY_FLAG, &
          nrow, mcol, nrow89,ncol89, ierr)
  ! TODO check all of these and find where they are called, nrow89,ncol89 perhaps in the invdist script the input to invdist should be modified
  
  ! Extra validation before proceeding
  write(LDT_logunit,*) '[DEBUG] After get_amsr_l1r - dimensions:'
  write(LDT_logunit,*) '   nrow=', nrow, ', mcol=', mcol
  
  ! Check array dimensions before calling L1RTB2ARFS_INVDIS
  if (nrow <= 0 .or. mcol <= 0) then
    write(LDT_logunit,*) '[ERR] Invalid array dimensions: nrow=', nrow, ', mcol=', mcol
    rc = 1
    return
  endif
  
  ! Verify lat/lon array sizes match brightness temperature array sizes
  if (size(LAT_L1R,1) /= size(TB_10H,1) .or. size(LAT_L1R,2) /= size(TB_10H,2) .or. &
      size(LON_L1R,1) /= size(TB_10H,1) .or. size(LON_L1R,2) /= size(TB_10H,2)) then
    write(LDT_logunit,*) '[ERR] Dimension mismatch between lat/lon and brightness temperature arrays'
    write(LDT_logunit,*) '      LAT_L1R: ', size(LAT_L1R,1), 'x', size(LAT_L1R,2)
    write(LDT_logunit,*) '      TB_10H: ', size(TB_10H,1), 'x', size(TB_10H,2)
    rc = 1
    return
  endif
  
  if (ierr == 1) then
     if (nrow == 0 .and. mcol == 0) then
        write(LDT_logunit,*)'[ERR] Problem reading ', trim(AMSRFILE)
        rc = 1
        return
     else
        write(LDT_logunit,*)'[ERR] Unknown internal error!'
        write(LDT_logunit,*)'[ERR] Aborting...'
        call LDT_endrun()
     end if
  end if
  
    
  CALL L1RTB2ARFS_INVDIS(TIME_L1R, TB_10H, TB_10V, TB_18H, TB_18V, TB_23H, TB_23V, &
          TB_36H, TB_36V, TB_89H, TB_89V, LAND_WATER_FRAC, &
          SNOW, PRECIP, QUALITY_FLAG, &  
          LAT_L1R, LON_L1R, nrow, mcol, &
          ARFS_LAT, ARFS_LON, ARFS_TIME, ARFS_LAND_WATER_FRAC, &
          ARFS_TB_10H, ARFS_TB_10V, ARFS_TB_18H, ARFS_TB_18V, &
          ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, &
          ARFS_TB_89H, ARFS_TB_89V, ARFS_QUALITY_FLAG, & 
          ARFS_SAMPLE_V, ARFS_SAMPLE_H)
  
  AMSReOPL%ARFS_TB_10H = ARFS_TB_10H
  AMSReOPL%ARFS_TB_10V = ARFS_TB_10V
  AMSReOPL%ARFS_TB_18H = ARFS_TB_18H
  AMSReOPL%ARFS_TB_18V = ARFS_TB_18V
  AMSReOPL%ARFS_TB_23H = ARFS_TB_23H
  AMSReOPL%ARFS_TB_23V = ARFS_TB_23V
  AMSReOPL%ARFS_TB_36H = ARFS_TB_36H
  AMSReOPL%ARFS_TB_36V = ARFS_TB_36V
  AMSReOPL%ARFS_TB_89H = ARFS_TB_89H
  AMSReOPL%ARFS_TB_89V = ARFS_TB_89V
  AMSReOPL%ARFS_LAND_WATER_FRAC = ARFS_LAND_WATER_FRAC
  AMSReOPL%ARFS_QUALITY_FLAG = ARFS_QUALITY_FLAG


    variable_name(1)  = 'ARFS_TIME'
    variable_name(2)  = 'ARFS_TB_10H'
    variable_name(3)  = 'ARFS_TB_10V'
    variable_name(4)  = 'ARFS_TB_18H'
    variable_name(5)  = 'ARFS_TB_18V'
    variable_name(6)  = 'ARFS_TB_23H'
    variable_name(7)  = 'ARFS_TB_23V'
    variable_name(8)  = 'ARFS_TB_36H'
    variable_name(9)  = 'ARFS_TB_36V'
    variable_name(10) = 'ARFS_TB_89H'
    variable_name(11) = 'ARFS_TB_89V'
    variable_name(12) = 'ARFS_LAND_WATER_FRAC'
    variable_name(13) = 'ARFS_QUALITY_FLAG'
    
    L1R_dir_len = len_trim(L1R_dir)
    L1R_fname_len = len_trim(AMSRFILE)
    Orbit = trim(AMSRFILE(L1R_dir_len+24:L1R_dir_len+24)) !E.J:  based on AMSR half-orbit file naming convention

    !if(AMSReOPL%L1Rtype.eq.1) then  !NRT 
      !Orbit = trim(AMSRFILE(L1R_dir_len+24:L1R_dir_len+24)) !E.J: this should be modified based on AMSR half-orbit files
    !elseif(AMSReOPL%L1Rtype.eq.2) then  !Historical
     !Orbit = trim(AMSRFILE(L1R_dir_len+20:L1R_dir_len+20))
    !endif
  
    !=================================================
    ! TODO: check the SMAPL1BTOL1C_ARFS.F90 for writing the file and modify the following for AMSR2 L1R reader
    ! Modified section for AMSR_L1R_ARFS_RESAMPLE.F90
    ! Replace the existing file writing section with this NetCDF version

    if(AMSReOPL%L1RresampWriteOpt.eq.1) then

        filename_start_pos = index(AMSRFILE, '/', back=.true.) + 1
        basename = AMSRFILE(filename_start_pos:)
        extracted_part = basename(8:len_trim(basename)-3)
       
       ! Construct NetCDF filename
       if(AMSReOPL%L1Rtype.eq.1) then  !NRT
          netcdf_filename = trim(AMSReOPL%L1Rresampledir_02)//"/AMSR_L1R_resampled_"//&
                            trim(AMSRFILE(L1R_dir_len+18:L1R_fname_len-3))//".nc"
       elseif(AMSReOPL%L1Rtype.eq.2) then  !Historical  
          netcdf_filename = trim(AMSReOPL%L1Rresampledir_02)//'/AMSR_L1R_resampled_'//trim(extracted_part)//'.nc'
       endif
       
       ! Create NetCDF file
       call LDT_verify(nf90_create(trim(netcdf_filename), NF90_NETCDF4, ncid), &
            '[ERR] nf90_create failed for AMSR resampled output')
       
       ! Define dimensions
       call LDT_verify(nf90_def_dim(ncid, 'time', 1, time_dimid), &
            '[ERR] nf90_def_dim failed for time')
       call LDT_verify(nf90_def_dim(ncid, 'lat', arfs_nrow_lat, lat_dimid), &
            '[ERR] nf90_def_dim failed for lat')
       call LDT_verify(nf90_def_dim(ncid, 'lon', arfs_mcol_lon, lon_dimid), &
            '[ERR] nf90_def_dim failed for lon')
       
       ! Define coordinate variables
       call LDT_verify(nf90_def_var(ncid, 'time', NF90_DOUBLE, time_dimid, time_varid), &
            '[ERR] nf90_def_var failed for time')
       call LDT_verify(nf90_def_var(ncid, 'lat', NF90_FLOAT, lat_dimid, lat_varid), &
            '[ERR] nf90_def_var failed for lat')  
       call LDT_verify(nf90_def_var(ncid, 'lon', NF90_FLOAT, lon_dimid, lon_varid), &
            '[ERR] nf90_def_var failed for lon')
       
       ! Define data variables - all TB channels and land water fraction
       call LDT_verify(nf90_def_var(ncid, 'TB_10H', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_10h_varid), &
            '[ERR] nf90_def_var failed for TB_10H')
       call LDT_verify(nf90_def_var(ncid, 'TB_10V', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_10v_varid), &
            '[ERR] nf90_def_var failed for TB_10V')
       call LDT_verify(nf90_def_var(ncid, 'TB_18H', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_18h_varid), &
            '[ERR] nf90_def_var failed for TB_18H')
       call LDT_verify(nf90_def_var(ncid, 'TB_18V', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_18v_varid), &
            '[ERR] nf90_def_var failed for TB_18V')
       call LDT_verify(nf90_def_var(ncid, 'TB_23H', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_23h_varid), &
            '[ERR] nf90_def_var failed for TB_23H')
       call LDT_verify(nf90_def_var(ncid, 'TB_23V', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_23v_varid), &
            '[ERR] nf90_def_var failed for TB_23V')
       call LDT_verify(nf90_def_var(ncid, 'TB_36H', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_36h_varid), &
            '[ERR] nf90_def_var failed for TB_36H')
       call LDT_verify(nf90_def_var(ncid, 'TB_36V', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_36v_varid), &
            '[ERR] nf90_def_var failed for TB_36V')
       call LDT_verify(nf90_def_var(ncid, 'TB_89H', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_89h_varid), &
            '[ERR] nf90_def_var failed for TB_89H')
       call LDT_verify(nf90_def_var(ncid, 'TB_89V', NF90_FLOAT, &
            [lon_dimid, lat_dimid], tb_89v_varid), &
            '[ERR] nf90_def_var failed for TB_89V')
       call LDT_verify(nf90_def_var(ncid, 'LAND_WATER_FRAC', NF90_FLOAT, &
            [lon_dimid, lat_dimid], lwf_varid), &
            '[ERR] nf90_def_var failed for LAND_WATER_FRAC')
       call LDT_verify(nf90_def_var(ncid, 'QUALITY_FLAG', NF90_BYTE, &
            [lon_dimid, lat_dimid], qf_varid), &
            '[ERR] nf90_def_var failed for QUALITY_FLAG')
       
        ! Add variable attributes for coordinate variables
       call LDT_verify(nf90_put_att(ncid, time_varid, 'units', 'seconds since 1970-01-01T00:00:00Z'), &
            '[ERR] nf90_put_att failed for time units')
       call LDT_verify(nf90_put_att(ncid, time_varid, 'standard_name', 'time'), &
            '[ERR] nf90_put_att failed for time standard_name')
       call LDT_verify(nf90_put_att(ncid, time_varid, 'calendar', 'standard'), &
            '[ERR] nf90_put_att failed for time calendar')
       call LDT_verify(nf90_put_att(ncid, lat_varid, 'units', 'degrees_north'), &
            '[ERR] nf90_put_att failed for lat units')
       call LDT_verify(nf90_put_att(ncid, lat_varid, 'standard_name', 'latitude'), &
            '[ERR] nf90_put_att failed for lat standard_name')
       call LDT_verify(nf90_put_att(ncid, lon_varid, 'units', 'degrees_east'), &
            '[ERR] nf90_put_att failed for lon units')
       call LDT_verify(nf90_put_att(ncid, lon_varid, 'standard_name', 'longitude'), &
            '[ERR] nf90_put_att failed for lon standard_name')
       
       ! TB_10H attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_10h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_10H units')
       call LDT_verify(nf90_put_att(ncid, tb_10h_varid, 'long_name', 'Brightness Temperature 10.65 GHz H-pol'), &
            '[ERR] nf90_put_att failed for TB_10H long_name')
       call LDT_verify(nf90_put_att(ncid, tb_10h_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_10H _FillValue')
       
       ! TB_10V attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_10v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_10V units')
       call LDT_verify(nf90_put_att(ncid, tb_10v_varid, 'long_name', 'Brightness Temperature 10.65 GHz V-pol'), &
            '[ERR] nf90_put_att failed for TB_10V long_name')
       call LDT_verify(nf90_put_att(ncid, tb_10v_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_10V _FillValue')
       
       ! TB_18H attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_18h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_18H units')
       call LDT_verify(nf90_put_att(ncid, tb_18h_varid, 'long_name', 'Brightness Temperature 18.7 GHz H-pol'), &
            '[ERR] nf90_put_att failed for TB_18H long_name')
       call LDT_verify(nf90_put_att(ncid, tb_18h_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_18H _FillValue')
       
       ! TB_18V attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_18v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_18V units')
       call LDT_verify(nf90_put_att(ncid, tb_18v_varid, 'long_name', 'Brightness Temperature 18.7 GHz V-pol'), &
            '[ERR] nf90_put_att failed for TB_18V long_name')
       call LDT_verify(nf90_put_att(ncid, tb_18v_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_18V _FillValue')
       
       ! TB_23H attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_23h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_23H units')
       call LDT_verify(nf90_put_att(ncid, tb_23h_varid, 'long_name', 'Brightness Temperature 23.8 GHz H-pol'), &
            '[ERR] nf90_put_att failed for TB_23H long_name')
       call LDT_verify(nf90_put_att(ncid, tb_23h_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_23H _FillValue')
       
       ! TB_23V attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_23v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_23V units')
       call LDT_verify(nf90_put_att(ncid, tb_23v_varid, 'long_name', 'Brightness Temperature 23.8 GHz V-pol'), &
            '[ERR] nf90_put_att failed for TB_23V long_name')
       call LDT_verify(nf90_put_att(ncid, tb_23v_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_23V _FillValue')
       
       ! TB_36H attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_36h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_36H units')
       call LDT_verify(nf90_put_att(ncid, tb_36h_varid, 'long_name', 'Brightness Temperature 36.5 GHz H-pol'), &
            '[ERR] nf90_put_att failed for TB_36H long_name')
       call LDT_verify(nf90_put_att(ncid, tb_36h_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_36H _FillValue')
       
       ! TB_36V attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_36v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_36V units')
       call LDT_verify(nf90_put_att(ncid, tb_36v_varid, 'long_name', 'Brightness Temperature 36.5 GHz V-pol'), &
            '[ERR] nf90_put_att failed for TB_36V long_name')
       call LDT_verify(nf90_put_att(ncid, tb_36v_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_36V _FillValue')
       
       ! TB_89H attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_89h_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_89H units')
       call LDT_verify(nf90_put_att(ncid, tb_89h_varid, 'long_name', 'Brightness Temperature 89.0 GHz H-pol'), &
            '[ERR] nf90_put_att failed for TB_89H long_name')
       call LDT_verify(nf90_put_att(ncid, tb_89h_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_89H _FillValue')
       
       ! TB_89V attributes with fill value
       call LDT_verify(nf90_put_att(ncid, tb_89v_varid, 'units', 'K'), &
            '[ERR] nf90_put_att failed for TB_89V units')
       call LDT_verify(nf90_put_att(ncid, tb_89v_varid, 'long_name', 'Brightness Temperature 89.0 GHz V-pol'), &
            '[ERR] nf90_put_att failed for TB_89V long_name')
       call LDT_verify(nf90_put_att(ncid, tb_89v_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for TB_89V _FillValue')

       ! LAND_WATER_FRAC attributes with fill value
       call LDT_verify(nf90_put_att(ncid, lwf_varid, 'units', 'fraction'), &
            '[ERR] nf90_put_att failed for LAND_WATER_FRAC units')
       call LDT_verify(nf90_put_att(ncid, lwf_varid, 'long_name', 'Land Water Fraction'), &
            '[ERR] nf90_put_att failed for LAND_WATER_FRAC long_name')
       call LDT_verify(nf90_put_att(ncid, lwf_varid, '_FillValue', -9999.0), &
            '[ERR] nf90_put_att failed for LAND_WATER_FRAC _FillValue')

       ! QUALITY_FLAG attributes
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'units', 'dimensionless'), &
            '[ERR] nf90_put_att failed for QUALITY_FLAG units')
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'long_name', &
            'Quality flags: bit0=ocean, bit1=precip, bit2=snow'), &
            '[ERR] nf90_put_att failed for QUALITY_FLAG long_name')
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'flag_meanings', &
            'ocean precipitation snow'), &
            '[ERR] nf90_put_att failed for QUALITY_FLAG flag_meanings')
       call LDT_verify(nf90_put_att(ncid, qf_varid, 'flag_masks', [1, 2, 4]), &
            '[ERR] nf90_put_att failed for QUALITY_FLAG flag_masks')
       call LDT_verify(nf90_put_att(ncid, qf_varid, '_FillValue', int(-128, kind=1)), &
            '[ERR] nf90_put_att failed for QUALITY_FLAG _FillValue')
            
       ! Add global attributes
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.10'), &
            '[ERR] nf90_put_att failed for Conventions')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'title', &
            'AMSR L1R Resampled to ARFS Grid'), &
            '[ERR] nf90_put_att failed for title')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'institution', &
            'NASA GSFC Hydrological Sciences Laboratory'), &
            '[ERR] nf90_put_att failed for institution')
       call LDT_verify(nf90_put_att(ncid, NF90_GLOBAL, 'source_file', &
            trim(AMSRFILE)), '[ERR] nf90_put_att failed for source_file')
       
       ! End definition mode
       call LDT_verify(nf90_enddef(ncid), '[ERR] nf90_enddef failed')
       
       ! Write coordinate data  
       allocate(lats(arfs_nrow_lat))
       allocate(lons(arfs_mcol_lon))

       lats = real(ARFS_LAT)  ! cast from REAL*8 to default REAL for NF90_FLOAT
       lons = real(ARFS_LON)

       call LDT_verify(nf90_put_var(ncid, lat_varid, lats), &
            '[ERR] nf90_put_var failed for lats')
       call LDT_verify(nf90_put_var(ncid, lon_varid, lons), &
            '[ERR] nf90_put_var failed for lons')
             
       ! Find where filename starts (after the last slash)
       filename_start_pos = index(AMSRFILE, '/', back=.true.) + 1
       
       ! Extract datetime from filename: YYYYMMDDHHMM at positions 8-19 of filename
       ! Example: GW1AM2_202401011913_135D... -> 202401011913
       datetime_str = basename(8:19)
       
       ! Parse date/time components from YYYYMMDDHHMM
       read(datetime_str(1:4), '(I4)') year
       read(datetime_str(5:6), '(I2)') month  
       read(datetime_str(7:8), '(I2)') day
       read(datetime_str(9:10), '(I2)') hour
       read(datetime_str(11:12), '(I2)') minute
       
       ! Create ESMF time objects
       call ESMF_TimeSet(reference_time, yy=1970, mm=1, dd=1, h=0, m=0, s=0, rc=rc_time)
       call ESMF_TimeSet(file_time, yy=year, mm=month, dd=day, h=hour, m=minute, s=0, rc=rc_time)
       
       ! Calculate time difference in seconds since epoch
       time_diff = file_time - reference_time
       call ESMF_TimeIntervalGet(time_diff, s_r8=time_seconds, rc=rc_time)
       
       call LDT_verify(nf90_put_var(ncid, time_varid, time_seconds), &
            '[ERR] nf90_put_var failed for time')
       
       ! Write all brightness temperature and fraction data
       call LDT_verify(nf90_put_var(ncid, tb_10h_varid, ARFS_TB_10H), &
            '[ERR] nf90_put_var failed for TB_10H')
       call LDT_verify(nf90_put_var(ncid, tb_10v_varid, ARFS_TB_10V), &
            '[ERR] nf90_put_var failed for TB_10V')
       call LDT_verify(nf90_put_var(ncid, tb_18h_varid, ARFS_TB_18H), &
            '[ERR] nf90_put_var failed for TB_18H')
       call LDT_verify(nf90_put_var(ncid, tb_18v_varid, ARFS_TB_18V), &
            '[ERR] nf90_put_var failed for TB_18V')
       call LDT_verify(nf90_put_var(ncid, tb_23h_varid, ARFS_TB_23H), &
            '[ERR] nf90_put_var failed for TB_23H')
       call LDT_verify(nf90_put_var(ncid, tb_23v_varid, ARFS_TB_23V), &
            '[ERR] nf90_put_var failed for TB_23V')
       call LDT_verify(nf90_put_var(ncid, tb_36h_varid, ARFS_TB_36H), &
            '[ERR] nf90_put_var failed for TB_36H')
       call LDT_verify(nf90_put_var(ncid, tb_36v_varid, ARFS_TB_36V), &
            '[ERR] nf90_put_var failed for TB_36V')
       call LDT_verify(nf90_put_var(ncid, tb_89h_varid, ARFS_TB_89H), &
            '[ERR] nf90_put_var failed for TB_89H')
       call LDT_verify(nf90_put_var(ncid, tb_89v_varid, ARFS_TB_89V), &
            '[ERR] nf90_put_var failed for TB_89V')
       call LDT_verify(nf90_put_var(ncid, lwf_varid, ARFS_LAND_WATER_FRAC), &
            '[ERR] nf90_put_var failed for LAND_WATER_FRAC')
       call LDT_verify(nf90_put_var(ncid, qf_varid, ARFS_QUALITY_FLAG), &
            '[ERR] nf90_put_var failed for QUALITY_FLAG')
            
       ! Close the file
       call LDT_verify(nf90_close(ncid), '[ERR] nf90_close failed')
       
       deallocate(lats)
       deallocate(lons)
       
       write(LDT_logunit,*) '[INFO] Successfully wrote NetCDF resampled file: ', trim(netcdf_filename)
       
    endif

    ! end of TODO for writting the outputfile
    !=================================================

  ! Cleanup [EJ: should I delocate the ARFS files?]
  deallocate(TB_10H)
  deallocate(TB_10V)
  deallocate(TB_18H)
  deallocate(TB_18V)
  deallocate(TB_23H)
  deallocate(TB_23V)
  deallocate(TB_36H)
  deallocate(TB_36V)
  deallocate(TB_89H)
  deallocate(TB_89V)
  
  deallocate(LAT_L1R)
  deallocate(LON_L1R)
  deallocate(LAND_WATER_FRAC)
  deallocate(RFI_FLAG)
  deallocate(SNOW)
  deallocate(PRECIP)
  
  deallocate(ARFS_LAT)
  deallocate(ARFS_LON)
  deallocate(QUALITY_FLAG)
end subroutine AMSR_L1R_RESAMPLE
