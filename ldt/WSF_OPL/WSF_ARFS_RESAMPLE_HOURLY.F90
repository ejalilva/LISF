!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! SUBROUTINE: WSF_ARFS_RESAMPLE_HOURLY
!
! DESCRIPTION: Simplified hourly resampling with proper band quality handling
!
!-------------------------------------------------------------------------

subroutine WSF_ARFS_RESAMPLE_HOURLY(hour_files, n_files, output_dir, &
                                    yyyymmdd, hour_str, n)

    USE TOOLSUBS_WSF
    USE invdist_wsf2arfs
    USE LDT_logMod
    USE LDT_coreMod, only: LDT_rc
    USE LDT_WSF_ARFS_netcdfMod, only: LDT_WSF_ARFS_write_netcdf_hourly
    USE LDT_wsf_oplMod, only: wsf_file_info

    IMPLICIT NONE
    
    ! Arguments
    integer, intent(in) :: n_files
    type(wsf_file_info), intent(in) :: hour_files(n_files)
    character(len=*), intent(in) :: output_dir
    character(len=8), intent(in) :: yyyymmdd
    character(len=2), intent(in) :: hour_str
    integer, intent(in) :: n
    
    ! ARFS grid arrays
    real*8, allocatable :: ARFS_LAT(:), ARFS_LON(:)
    
    ! Accumulated data arrays
    real*4, allocatable :: ARFS_TB_10H_SUM(:,:), ARFS_TB_10V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_18H_SUM(:,:), ARFS_TB_18V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_23H_SUM(:,:), ARFS_TB_23V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_36H_SUM(:,:), ARFS_TB_36V_SUM(:,:)
    real*4, allocatable :: ARFS_TB_89H_SUM(:,:), ARFS_TB_89V_SUM(:,:)
    real*4, allocatable :: ARFS_LAND_FRAC_SUM(:,:)
    
    ! Count arrays
    integer*4, allocatable :: ARFS_COUNT_10H(:,:), ARFS_COUNT_10V(:,:)
    integer*4, allocatable :: ARFS_COUNT_18H(:,:), ARFS_COUNT_18V(:,:)
    integer*4, allocatable :: ARFS_COUNT_23H(:,:), ARFS_COUNT_23V(:,:)
    integer*4, allocatable :: ARFS_COUNT_36H(:,:), ARFS_COUNT_36V(:,:)
    integer*4, allocatable :: ARFS_COUNT_89H(:,:), ARFS_COUNT_89V(:,:)
    integer*4, allocatable :: ARFS_COUNT_LAND(:,:)
    
    ! Final arrays
    real*8, allocatable :: ARFS_TIME(:,:)
    real*4, allocatable :: ARFS_TB_10H(:,:), ARFS_TB_10V(:,:)
    real*4, allocatable :: ARFS_TB_18H(:,:), ARFS_TB_18V(:,:)
    real*4, allocatable :: ARFS_TB_23H(:,:), ARFS_TB_23V(:,:)
    real*4, allocatable :: ARFS_TB_36H(:,:), ARFS_TB_36V(:,:)
    real*4, allocatable :: ARFS_TB_89H(:,:), ARFS_TB_89V(:,:)
    real*4, allocatable :: ARFS_LAND_FRAC(:,:)
    integer*1, allocatable :: ARFS_QUALITY_FLAG(:,:)
    integer*4, allocatable :: ARFS_SAMPLE_V(:,:), ARFS_SAMPLE_H(:,:)
    
    ! Temporary arrays for each file
    real*8, allocatable :: TEMP_TIME(:,:)
    real*4, allocatable :: TEMP_TB_10H(:,:), TEMP_TB_10V(:,:)
    real*4, allocatable :: TEMP_TB_18H(:,:), TEMP_TB_18V(:,:)
    real*4, allocatable :: TEMP_TB_23H(:,:), TEMP_TB_23V(:,:)
    real*4, allocatable :: TEMP_TB_36H(:,:), TEMP_TB_36V(:,:)
    real*4, allocatable :: TEMP_TB_89H(:,:), TEMP_TB_89V(:,:)
    real*4, allocatable :: TEMP_LAND_FRAC(:,:)
    integer*1, allocatable :: TEMP_QUALITY_FLAG(:,:)
    integer*4, allocatable :: TEMP_SAMPLE_V(:,:), TEMP_SAMPLE_H(:,:)
    
    ! WSF input data
    real*4, allocatable :: tb_lowres(:,:,:)
    real*4, allocatable :: lat_in(:,:)
    real*4, allocatable :: lon_in(:,:)
    real*4, allocatable :: land_frac_low(:,:)
    integer*4, allocatable :: quality_flag_in(:,:)
    integer*4, allocatable :: band_quality_flags(:,:,:)  ! NEW
    real*4, allocatable :: earth_inc_angle(:,:,:)
    integer*4, allocatable :: snow_in(:,:)
    integer*4, allocatable :: precip_in(:,:)
    
    ! Channel extraction arrays
    real*4, allocatable :: tb_10h(:,:), tb_10v(:,:)
    real*4, allocatable :: tb_18h(:,:), tb_18v(:,:)
    real*4, allocatable :: tb_23h(:,:), tb_23v(:,:)
    real*4, allocatable :: tb_36h(:,:), tb_36v(:,:)
    real*4, allocatable :: tb_89h(:,:), tb_89v(:,:)
    real*8, allocatable :: time_array(:)
    
    real*4, allocatable :: chan_frequencies(:)
    character*1, allocatable :: chan_polarizations(:)
    
    ! Quality flag accumulation
    integer*4, allocatable :: ARFS_QF_COUNTS(:,:,:)  ! (2560, 1920, 8 bits)
    integer*4, allocatable :: ARFS_QF_TOTALS(:,:)
    
    integer :: i, j, r, c, ifile
    integer :: nscans, nfovs, nchans, ierr, ichan
    real :: freq
    character*1 :: pol
    character(len=255) :: output_filename
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] WSF HOURLY GROUP PROCESSING - FIXED'
    write(LDT_logunit,*)'[INFO] Hour: ', hour_str, 'H'
    write(LDT_logunit,*)'[INFO] Number of files: ', n_files
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Setup ARFS grid
    allocate(ARFS_LAT(1920))
    allocate(ARFS_LON(2560))
    
    do r = 1, 1920
        ARFS_LAT(r) = 90.0 - (r - 0.5) * 0.09375
    end do
    
    do c = 1, 2560
        ARFS_LON(c) = -180.0 + (c - 0.5) * 0.140625
    end do
    
    ! Allocate accumulation arrays
    allocate(ARFS_TB_10H_SUM(2560,1920), ARFS_TB_10V_SUM(2560,1920))
    allocate(ARFS_TB_18H_SUM(2560,1920), ARFS_TB_18V_SUM(2560,1920))
    allocate(ARFS_TB_23H_SUM(2560,1920), ARFS_TB_23V_SUM(2560,1920))
    allocate(ARFS_TB_36H_SUM(2560,1920), ARFS_TB_36V_SUM(2560,1920))
    allocate(ARFS_TB_89H_SUM(2560,1920), ARFS_TB_89V_SUM(2560,1920))
    allocate(ARFS_LAND_FRAC_SUM(2560,1920))
    
    allocate(ARFS_COUNT_10H(2560,1920), ARFS_COUNT_10V(2560,1920))
    allocate(ARFS_COUNT_18H(2560,1920), ARFS_COUNT_18V(2560,1920))
    allocate(ARFS_COUNT_23H(2560,1920), ARFS_COUNT_23V(2560,1920))
    allocate(ARFS_COUNT_36H(2560,1920), ARFS_COUNT_36V(2560,1920))
    allocate(ARFS_COUNT_89H(2560,1920), ARFS_COUNT_89V(2560,1920))
    allocate(ARFS_COUNT_LAND(2560,1920))
    
    allocate(TEMP_TIME(2560,1920))
    allocate(TEMP_TB_10H(2560,1920), TEMP_TB_10V(2560,1920))
    allocate(TEMP_TB_18H(2560,1920), TEMP_TB_18V(2560,1920))
    allocate(TEMP_TB_23H(2560,1920), TEMP_TB_23V(2560,1920))
    allocate(TEMP_TB_36H(2560,1920), TEMP_TB_36V(2560,1920))
    allocate(TEMP_TB_89H(2560,1920), TEMP_TB_89V(2560,1920))
    allocate(TEMP_LAND_FRAC(2560,1920))
    allocate(TEMP_QUALITY_FLAG(2560,1920))
    allocate(TEMP_SAMPLE_V(2560,1920), TEMP_SAMPLE_H(2560,1920))
    
    allocate(ARFS_QF_COUNTS(2560,1920,8))
    allocate(ARFS_QF_TOTALS(2560,1920))
    
    ! Initialize accumulation arrays
    ARFS_TB_10H_SUM = 0.0
    ARFS_TB_10V_SUM = 0.0
    ARFS_TB_18H_SUM = 0.0
    ARFS_TB_18V_SUM = 0.0
    ARFS_TB_23H_SUM = 0.0
    ARFS_TB_23V_SUM = 0.0
    ARFS_TB_36H_SUM = 0.0
    ARFS_TB_36V_SUM = 0.0
    ARFS_TB_89H_SUM = 0.0
    ARFS_TB_89V_SUM = 0.0
    ARFS_LAND_FRAC_SUM = 0.0
    
    ARFS_COUNT_10H = 0
    ARFS_COUNT_10V = 0
    ARFS_COUNT_18H = 0
    ARFS_COUNT_18V = 0
    ARFS_COUNT_23H = 0
    ARFS_COUNT_23V = 0
    ARFS_COUNT_36H = 0
    ARFS_COUNT_36V = 0
    ARFS_COUNT_89H = 0
    ARFS_COUNT_89V = 0
    ARFS_COUNT_LAND = 0
    ARFS_QF_COUNTS = 0
    ARFS_QF_TOTALS = 0
    
    ! =====================================================================
    ! PROCESS EACH FILE
    ! =====================================================================
    do ifile = 1, n_files
        write(LDT_logunit,*)'[INFO] ----------------------------------------'
        write(LDT_logunit,*)'[INFO] Processing file ', ifile, '/', n_files
        write(LDT_logunit,*)'[INFO] ', trim(hour_files(ifile)%filename)
        
        ! Initialize temporary arrays
        TEMP_TIME = 0.0
        TEMP_TB_10H = 0.0
        TEMP_TB_10V = 0.0
        TEMP_TB_18H = 0.0
        TEMP_TB_18V = 0.0
        TEMP_TB_23H = 0.0
        TEMP_TB_23V = 0.0
        TEMP_TB_36H = 0.0
        TEMP_TB_36V = 0.0
        TEMP_TB_89H = 0.0
        TEMP_TB_89V = 0.0
        TEMP_LAND_FRAC = 0.0
        TEMP_QUALITY_FLAG = 0
        TEMP_SAMPLE_V = 0
        TEMP_SAMPLE_H = 0
        
        ! Read WSF data WITH band quality flags
        call get_wsf_data_with_flags(hour_files(ifile)%filename, &
            tb_lowres, lat_in, lon_in, land_frac_low, quality_flag_in, &
            band_quality_flags, &  ! NEW parameter
            earth_inc_angle, snow_in, precip_in, &
            nscans, nfovs, nchans, &
            chan_frequencies, chan_polarizations, &
            ierr)
        
        if (ierr /= 0) then
            write(LDT_logunit,*)'[WARN] Failed to read file, skipping'
            cycle
        endif
        
        write(LDT_logunit,*)'[INFO] Read data dimensions:'
        write(LDT_logunit,*)'[INFO]   nFOVR  = ', nfovs
        write(LDT_logunit,*)'[INFO]   nScanR = ', nscans
        write(LDT_logunit,*)'[INFO]   nChan  = ', nchans
        
        ! Extract channels - transpose from (nfovs, nscans) to (nscans, nfovs)
        allocate(tb_10h(nscans, nfovs))
        allocate(tb_10v(nscans, nfovs))
        allocate(tb_18h(nscans, nfovs))
        allocate(tb_18v(nscans, nfovs))
        allocate(tb_23h(nscans, nfovs))
        allocate(tb_23v(nscans, nfovs))
        allocate(tb_36h(nscans, nfovs))
        allocate(tb_36v(nscans, nfovs))
        allocate(tb_89h(nscans, nfovs))
        allocate(tb_89v(nscans, nfovs))
        allocate(time_array(nscans))
        
        ! Initialize to fill values
        tb_10h = -9999.0
        tb_10v = -9999.0
        tb_18h = -9999.0
        tb_18v = -9999.0
        tb_23h = -9999.0
        tb_23v = -9999.0
        tb_36h = -9999.0
        tb_36v = -9999.0
        tb_89h = -9999.0
        tb_89v = -9999.0
        
        ! Extract each channel
        do ichan = 1, nchans
            freq = chan_frequencies(ichan)
            pol = chan_polarizations(ichan)
            
            ! Transpose while extracting
            if (abs(freq - 10.65) < 1.0) then
                if (pol == 'v' .or. pol == 'V') then
                    tb_10v = TRANSPOSE(tb_lowres(:,:,ichan))
                else if (pol == 'h' .or. pol == 'H') then
                    tb_10h = TRANSPOSE(tb_lowres(:,:,ichan))
                endif
            else if (abs(freq - 18.7) < 1.0) then
                if (pol == 'v' .or. pol == 'V') then
                    tb_18v = TRANSPOSE(tb_lowres(:,:,ichan))
                else if (pol == 'h' .or. pol == 'H') then
                    tb_18h = TRANSPOSE(tb_lowres(:,:,ichan))
                endif
            else if (abs(freq - 23.8) < 1.0) then
                if (pol == 'v' .or. pol == 'V') then
                    tb_23v = TRANSPOSE(tb_lowres(:,:,ichan))
                else if (pol == 'h' .or. pol == 'H') then
                    tb_23h = TRANSPOSE(tb_lowres(:,:,ichan))
                endif
            else if (abs(freq - 36.5) < 1.0) then
                if (pol == 'v' .or. pol == 'V') then
                    tb_36v = TRANSPOSE(tb_lowres(:,:,ichan))
                else if (pol == 'h' .or. pol == 'H') then
                    tb_36h = TRANSPOSE(tb_lowres(:,:,ichan))
                endif
            else if (abs(freq - 89.0) < 1.0) then
                if (pol == 'v' .or. pol == 'V') then
                    tb_89v = TRANSPOSE(tb_lowres(:,:,ichan))
                else if (pol == 'h' .or. pol == 'H') then
                    tb_89h = TRANSPOSE(tb_lowres(:,:,ichan))
                endif
            endif
        end do
        
        time_array = 0.0
        
        write(LDT_logunit,*)'[INFO] Calling inverse distance resampler WITH band quality...'
        
        ! Call resampler WITH band quality flags
        call WSF2ARFS_INVDIS(time_array, &
            tb_10h, tb_10v, tb_18h, tb_18v, &
            tb_23h, tb_23v, tb_36h, tb_36v, &
            tb_89h, tb_89v, &
            TRANSPOSE(land_frac_low), &
            TRANSPOSE(snow_in), &
            TRANSPOSE(precip_in), &
            TRANSPOSE(quality_flag_in), &
            reshape(band_quality_flags, (/nscans, nfovs, 6/), order=(/2,1,3/)), &  ! NEW: transpose band flags
            TRANSPOSE(lat_in), &
            TRANSPOSE(lon_in), &
            nscans, nfovs, &
            ARFS_LAT, ARFS_LON, &
            TEMP_TIME, TEMP_LAND_FRAC, &
            TEMP_TB_10H, TEMP_TB_10V, &
            TEMP_TB_18H, TEMP_TB_18V, &
            TEMP_TB_23H, TEMP_TB_23V, &
            TEMP_TB_36H, TEMP_TB_36V, &
            TEMP_TB_89H, TEMP_TB_89V, &
            TEMP_QUALITY_FLAG, &
            TEMP_SAMPLE_V, TEMP_SAMPLE_H)
        
        ! Accumulate data
        do r = 1, 1920
            do c = 1, 2560
                ! TB accumulation
                if (TEMP_TB_10H(c,r) > 0.0) then
                    ARFS_TB_10H_SUM(c,r) = ARFS_TB_10H_SUM(c,r) + TEMP_TB_10H(c,r)
                    ARFS_COUNT_10H(c,r) = ARFS_COUNT_10H(c,r) + 1
                endif
                if (TEMP_TB_10V(c,r) > 0.0) then
                    ARFS_TB_10V_SUM(c,r) = ARFS_TB_10V_SUM(c,r) + TEMP_TB_10V(c,r)
                    ARFS_COUNT_10V(c,r) = ARFS_COUNT_10V(c,r) + 1
                endif
                if (TEMP_TB_18H(c,r) > 0.0) then
                    ARFS_TB_18H_SUM(c,r) = ARFS_TB_18H_SUM(c,r) + TEMP_TB_18H(c,r)
                    ARFS_COUNT_18H(c,r) = ARFS_COUNT_18H(c,r) + 1
                endif
                if (TEMP_TB_18V(c,r) > 0.0) then
                    ARFS_TB_18V_SUM(c,r) = ARFS_TB_18V_SUM(c,r) + TEMP_TB_18V(c,r)
                    ARFS_COUNT_18V(c,r) = ARFS_COUNT_18V(c,r) + 1
                endif
                if (TEMP_TB_23H(c,r) > 0.0) then
                    ARFS_TB_23H_SUM(c,r) = ARFS_TB_23H_SUM(c,r) + TEMP_TB_23H(c,r)
                    ARFS_COUNT_23H(c,r) = ARFS_COUNT_23H(c,r) + 1
                endif
                if (TEMP_TB_23V(c,r) > 0.0) then
                    ARFS_TB_23V_SUM(c,r) = ARFS_TB_23V_SUM(c,r) + TEMP_TB_23V(c,r)
                    ARFS_COUNT_23V(c,r) = ARFS_COUNT_23V(c,r) + 1
                endif
                if (TEMP_TB_36H(c,r) > 0.0) then
                    ARFS_TB_36H_SUM(c,r) = ARFS_TB_36H_SUM(c,r) + TEMP_TB_36H(c,r)
                    ARFS_COUNT_36H(c,r) = ARFS_COUNT_36H(c,r) + 1
                endif
                if (TEMP_TB_36V(c,r) > 0.0) then
                    ARFS_TB_36V_SUM(c,r) = ARFS_TB_36V_SUM(c,r) + TEMP_TB_36V(c,r)
                    ARFS_COUNT_36V(c,r) = ARFS_COUNT_36V(c,r) + 1
                endif
                if (TEMP_TB_89H(c,r) > 0.0) then
                    ARFS_TB_89H_SUM(c,r) = ARFS_TB_89H_SUM(c,r) + TEMP_TB_89H(c,r)
                    ARFS_COUNT_89H(c,r) = ARFS_COUNT_89H(c,r) + 1
                endif
                if (TEMP_TB_89V(c,r) > 0.0) then
                    ARFS_TB_89V_SUM(c,r) = ARFS_TB_89V_SUM(c,r) + TEMP_TB_89V(c,r)
                    ARFS_COUNT_89V(c,r) = ARFS_COUNT_89V(c,r) + 1
                endif
                if (TEMP_LAND_FRAC(c,r) > 0.0) then
                    ARFS_LAND_FRAC_SUM(c,r) = ARFS_LAND_FRAC_SUM(c,r) + TEMP_LAND_FRAC(c,r)
                    ARFS_COUNT_LAND(c,r) = ARFS_COUNT_LAND(c,r) + 1
                endif
                
                ! Accumulate quality flags
                if (INT(TEMP_QUALITY_FLAG(c,r)) >= 0) then
                    ARFS_QF_TOTALS(c,r) = ARFS_QF_TOTALS(c,r) + 1
                    do i = 0, 7
                        if (IBITS(TEMP_QUALITY_FLAG(c,r), i, 1) == 1) then
                            ARFS_QF_COUNTS(c,r,i+1) = ARFS_QF_COUNTS(c,r,i+1) + 1
                        endif
                    end do
                endif
                
                ! Accumulate sample counts
                ARFS_SAMPLE_V(c,r) = ARFS_SAMPLE_V(c,r) + TEMP_SAMPLE_V(c,r)
                ARFS_SAMPLE_H(c,r) = ARFS_SAMPLE_H(c,r) + TEMP_SAMPLE_H(c,r)
            end do
        end do
        
        ! Cleanup for this file
        deallocate(tb_lowres, lat_in, lon_in)
        deallocate(land_frac_low, quality_flag_in, band_quality_flags)
        deallocate(earth_inc_angle, snow_in, precip_in)
        deallocate(chan_frequencies, chan_polarizations)
        deallocate(tb_10h, tb_10v, tb_18h, tb_18v)
        deallocate(tb_23h, tb_23v, tb_36h, tb_36v)
        deallocate(tb_89h, tb_89v, time_array)
        
    end do ! ifile loop
    
    ! =====================================================================
    ! CALCULATE MEANS
    ! =====================================================================
    allocate(ARFS_TIME(2560,1920))
    allocate(ARFS_TB_10H(2560,1920), ARFS_TB_10V(2560,1920))
    allocate(ARFS_TB_18H(2560,1920), ARFS_TB_18V(2560,1920))
    allocate(ARFS_TB_23H(2560,1920), ARFS_TB_23V(2560,1920))
    allocate(ARFS_TB_36H(2560,1920), ARFS_TB_36V(2560,1920))
    allocate(ARFS_TB_89H(2560,1920), ARFS_TB_89V(2560,1920))
    allocate(ARFS_LAND_FRAC(2560,1920))
    allocate(ARFS_QUALITY_FLAG(2560,1920))
    
    ARFS_TIME = 0.0
    
    ! Calculate means
    do r = 1, 1920
        do c = 1, 2560
            if (ARFS_COUNT_10H(c,r) > 0) then
                ARFS_TB_10H(c,r) = ARFS_TB_10H_SUM(c,r) / real(ARFS_COUNT_10H(c,r))
            else
                ARFS_TB_10H(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_10V(c,r) > 0) then
                ARFS_TB_10V(c,r) = ARFS_TB_10V_SUM(c,r) / real(ARFS_COUNT_10V(c,r))
            else
                ARFS_TB_10V(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_18H(c,r) > 0) then
                ARFS_TB_18H(c,r) = ARFS_TB_18H_SUM(c,r) / real(ARFS_COUNT_18H(c,r))
            else
                ARFS_TB_18H(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_18V(c,r) > 0) then
                ARFS_TB_18V(c,r) = ARFS_TB_18V_SUM(c,r) / real(ARFS_COUNT_18V(c,r))
            else
                ARFS_TB_18V(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_23H(c,r) > 0) then
                ARFS_TB_23H(c,r) = ARFS_TB_23H_SUM(c,r) / real(ARFS_COUNT_23H(c,r))
            else
                ARFS_TB_23H(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_23V(c,r) > 0) then
                ARFS_TB_23V(c,r) = ARFS_TB_23V_SUM(c,r) / real(ARFS_COUNT_23V(c,r))
            else
                ARFS_TB_23V(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_36H(c,r) > 0) then
                ARFS_TB_36H(c,r) = ARFS_TB_36H_SUM(c,r) / real(ARFS_COUNT_36H(c,r))
            else
                ARFS_TB_36H(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_36V(c,r) > 0) then
                ARFS_TB_36V(c,r) = ARFS_TB_36V_SUM(c,r) / real(ARFS_COUNT_36V(c,r))
            else
                ARFS_TB_36V(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_89H(c,r) > 0) then
                ARFS_TB_89H(c,r) = ARFS_TB_89H_SUM(c,r) / real(ARFS_COUNT_89H(c,r))
            else
                ARFS_TB_89H(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_89V(c,r) > 0) then
                ARFS_TB_89V(c,r) = ARFS_TB_89V_SUM(c,r) / real(ARFS_COUNT_89V(c,r))
            else
                ARFS_TB_89V(c,r) = -9999.0
            endif
            
            if (ARFS_COUNT_LAND(c,r) > 0) then
                ARFS_LAND_FRAC(c,r) = ARFS_LAND_FRAC_SUM(c,r) / real(ARFS_COUNT_LAND(c,r))
            else
                ARFS_LAND_FRAC(c,r) = -9999.0
            endif
            
            ! Quality flag via majority voting
            if (ARFS_QF_TOTALS(c,r) > 0) then
                ARFS_QUALITY_FLAG(c,r) = 0
                do i = 1, 8
                    if (ARFS_QF_COUNTS(c,r,i) > ARFS_QF_TOTALS(c,r)/2) then
                        ARFS_QUALITY_FLAG(c,r) = IOR(ARFS_QUALITY_FLAG(c,r), 2**(i-1))
                    endif
                end do
            else
                ARFS_QUALITY_FLAG(c,r) = -1
            endif
        end do
    end do
    
    ! =====================================================================
    ! WRITE OUTPUT
    ! =====================================================================
    output_filename = trim(output_dir)//'/WSF_SDR_resampled_'//yyyymmdd//'_t'//hour_str//'00.nc'
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Writing stitched hourly output'
    write(LDT_logunit,*)'[INFO] Output file: ', trim(output_filename)
    
    call LDT_WSF_ARFS_write_netcdf_hourly(2560, 1920, &
        ARFS_TB_10H, ARFS_TB_10V, ARFS_TB_18H, ARFS_TB_18V, &
        ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, &
        ARFS_TB_89H, ARFS_TB_89V, ARFS_LAND_FRAC, &
        ARFS_QUALITY_FLAG, ARFS_SAMPLE_V, ARFS_SAMPLE_H, &
        ARFS_LAT, ARFS_LON, output_filename, &
        yyyymmdd, hour_str, n_files)
    
    write(LDT_logunit,*)'[INFO] ========================================='
    write(LDT_logunit,*)'[INFO] Hourly stitching statistics:'
    write(LDT_logunit,*)'[INFO] Files processed: ', n_files
    write(LDT_logunit,*)'[INFO] Grid points with data:'
    write(LDT_logunit,*)'[INFO]   10H: ', count(ARFS_TB_10H > 0.0)
    write(LDT_logunit,*)'[INFO]   10V: ', count(ARFS_TB_10V > 0.0)
    write(LDT_logunit,*)'[INFO] Quality flag bits set:'
    write(LDT_logunit,*)'[INFO]   Ocean (bit 0): ', count(IBITS(ARFS_QUALITY_FLAG, 0, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Precip (bit 1): ', count(IBITS(ARFS_QUALITY_FLAG, 1, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Snow (bit 2): ', count(IBITS(ARFS_QUALITY_FLAG, 2, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 10GHz (bit 3): ', count(IBITS(ARFS_QUALITY_FLAG, 3, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 18GHz (bit 4): ', count(IBITS(ARFS_QUALITY_FLAG, 4, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 23GHz (bit 5): ', count(IBITS(ARFS_QUALITY_FLAG, 5, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 36GHz (bit 6): ', count(IBITS(ARFS_QUALITY_FLAG, 6, 1) == 1)
    write(LDT_logunit,*)'[INFO]   Bad 89GHz (bit 7): ', count(IBITS(ARFS_QUALITY_FLAG, 7, 1) == 1)
    write(LDT_logunit,*)'[INFO] ========================================='
    
    ! Cleanup
    deallocate(ARFS_LAT, ARFS_LON)
    deallocate(ARFS_TB_10H_SUM, ARFS_TB_10V_SUM)
    deallocate(ARFS_TB_18H_SUM, ARFS_TB_18V_SUM)
    deallocate(ARFS_TB_23H_SUM, ARFS_TB_23V_SUM)
    deallocate(ARFS_TB_36H_SUM, ARFS_TB_36V_SUM)
    deallocate(ARFS_TB_89H_SUM, ARFS_TB_89V_SUM)
    deallocate(ARFS_LAND_FRAC_SUM)
    deallocate(ARFS_COUNT_10H, ARFS_COUNT_10V)
    deallocate(ARFS_COUNT_18H, ARFS_COUNT_18V)
    deallocate(ARFS_COUNT_23H, ARFS_COUNT_23V)
    deallocate(ARFS_COUNT_36H, ARFS_COUNT_36V)
    deallocate(ARFS_COUNT_89H, ARFS_COUNT_89V)
    deallocate(ARFS_COUNT_LAND)
    deallocate(ARFS_QF_COUNTS, ARFS_QF_TOTALS)
    deallocate(TEMP_TIME, TEMP_TB_10H, TEMP_TB_10V)
    deallocate(TEMP_TB_18H, TEMP_TB_18V)
    deallocate(TEMP_TB_23H, TEMP_TB_23V)
    deallocate(TEMP_TB_36H, TEMP_TB_36V)
    deallocate(TEMP_TB_89H, TEMP_TB_89V)
    deallocate(TEMP_LAND_FRAC, TEMP_QUALITY_FLAG)
    deallocate(TEMP_SAMPLE_V, TEMP_SAMPLE_H)
    deallocate(ARFS_TIME, ARFS_TB_10H, ARFS_TB_10V)
    deallocate(ARFS_TB_18H, ARFS_TB_18V)
    deallocate(ARFS_TB_23H, ARFS_TB_23V)
    deallocate(ARFS_TB_36H, ARFS_TB_36V)
    deallocate(ARFS_TB_89H, ARFS_TB_89V)
    deallocate(ARFS_LAND_FRAC, ARFS_QUALITY_FLAG)

end subroutine WSF_ARFS_RESAMPLE_HOURLY
