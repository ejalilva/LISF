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
! MODULE: invdist_l1r2arfs
!
! REVISION HISTORY: 
!  26 Oct 2021: P.W.LIU; Initial implemetation
!  28 Jan 2022: P.W.LIU; CHNAGE OUTPUT COORDINATE TO COMPLY LIS OUTPUT
!  22 Feb 2022: Yonghwan Kwon; modified for LDT
!  10 Feb 2023: Eric Kemp; Added code to resample subset of variables.
!  10 Feb 2025: Ehsan Jalilvand; Modified the SMAP code to resmaple AMSR2 L1R
!
! DESCRIPTION: SUBROUTINE TO RESAMPLE L1R TB ONTO AIR FORCE GRID
!-------------------------------------------------------------------------

 MODULE invdist_l1r2arfs
   IMPLICIT NONE
 CONTAINS
   SUBROUTINE L1RTB2ARFS_INVDIS(tim, tb_10h, tb_10v, tb_18h, tb_18v, tb_23h, tb_23v, &
          tb_36h, tb_36v, tb_89h, tb_89v, land_water_frac, &
          snow_flag, precip_flag, quality_flag, & 
          lat_l1r, lon_l1r, nrows_l1rtb, ncols_l1rtb, &
          !lat89, lon89, nrows_89, ncols_89, &
          !rfi_flag, arfs_rfi_flag
          ref_lat, ref_lon, arfs_tim, arfs_land_water_frac, &
          arfs_tb_10h, arfs_tb_10v, arfs_tb_18h, arfs_tb_18v, &
          arfs_tb_23h, arfs_tb_23v, arfs_tb_36h, arfs_tb_36v, &
          arfs_tb_89h, arfs_tb_89v, arfs_quality_flag, &
          arfs_samplenumv, arfs_samplenumh)
   !SUBROUTINE L1BTB2ARFS_INVDIS(tim, tbvl1b_cor, tbhl1b_cor, tbvl1b, tbhl1b, surwat_v_l1b, surwat_h_l1b, &
        !netd_v_l1b, netd_h_l1b, lat_l1b, lon_l1b, tbv_qual_flag, tbh_qual_flag, sc_nadir_angle, antenna_scan_angle, nrows_l1btb, ncols_l1btb, &
        !ref_lat, ref_lon, arfs_tim, arfs_tbv_cor, arfs_tbh_cor, arfs_tbv, arfs_tbh, arfs_nedtv, arfs_nedth, &
        !arfs_surwatv, arfs_surwath, arfs_wt_cor_tbv, arfs_wt_cor_tbh, arfs_samplenumv, arfs_samplenumh)

     USE LDT_logMod, only: LDT_logunit ! Add this import for logging

     INTEGER(4) :: i, j, ii, jj, k, r, c, rr, rmin, rmax, cc, cmin, cmax, nrows_l1rtb, ncols_l1rtb
     INTEGER(4), PARAMETER :: qualitybit = 0
     REAL(8), PARAMETER :: RE_KM = 6371.228, search_radius = 20.0, PI = 3.141592653589793238, d2r = PI/180.0
     REAL(8)  :: gcdist, lat1, lon1, lat2, lon2
     LOGICAL :: has_snow, has_precip
     REAL*8,DIMENSION(nrows_l1rtb) :: tim
     REAL*4,DIMENSION(nrows_l1rtb,ncols_l1rtb) :: tb_10h, tb_10v, tb_18h, tb_18v, tb_23h, tb_23v, tb_36h, tb_36v, tb_89h, tb_89v
     REAL*4,DIMENSION(nrows_l1rtb,ncols_l1rtb) :: lat_l1r, lon_l1r
     INTEGER*4,DIMENSION(nrows_l1rtb,ncols_l1rtb) :: snow_flag, precip_flag, land_water_frac

     INTEGER*1,DIMENSION(nrows_l1rtb,ncols_l1rtb) :: quality_flag            ! Input from TOOLSUBS

     !! TODO: uncomment for RFI_flag implementation
     !INTEGER(4) :: nrows_89, ncols_89
     !REAL*4,DIMENSION(nrows_89,ncols_89) :: lat89, lon89
     !INTEGER*4,DIMENSION(nrows_89,ncols_89) :: rfi_flag
     !INTEGER*4,DIMENSION(2560,1920) :: arfs_rfi_flag
     !REAL*4,DIMENSION(2560,1920) :: arfs_wt_rfi_flag
     
     INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: zerodistflag
     REAL*8,DIMENSION(:),ALLOCATABLE :: ref_lat, ref_lon
     REAL*8,DIMENSION(2560,1920) :: arfs_tim
     REAL*4,DIMENSION(2560,1920) :: arfs_tb_10h, arfs_tb_10v, arfs_tb_18h, arfs_tb_18v, arfs_tb_23h, arfs_tb_23v, arfs_tb_36h, arfs_tb_36v, arfs_tb_89h, arfs_tb_89v, arfs_land_water_frac
     REAL*4,DIMENSION(2560,1920) :: arfs_wt_tim, arfs_wt_tb10v, arfs_wt_tb10h, arfs_wt_tb18v, arfs_wt_tb18h, arfs_wt_tb23v, arfs_wt_tb23h, arfs_wt_tb36v, arfs_wt_tb36h, arfs_wt_tb89v, arfs_wt_tb89h, arfs_wt_land_water_frac
     INTEGER*4,DIMENSION(2560,1920) :: arfs_samplenumv, arfs_samplenumh
     
     INTEGER*1,DIMENSION(2560,1920) :: arfs_quality_flag                     ! Output grid
     INTEGER,DIMENSION(2560,1920) :: snow_count, precip_count, ocean_count, total_count  ! Counters for majority vote
     INTEGER,DIMENSION(2560,1920) :: excluded_snow_count, excluded_precip_count


     !ALLOCATE(zerodistflag(size(ref_lat),size(ref_lon)))
     ALLOCATE(zerodistflag(size(ref_lon),size(ref_lat)))
     zerodistflag = 0 ! E.J
     
     !INITIAL THE OUTPUT VARIABLES WITH FILL VALUE
    arfs_tim=0.0
    arfs_tb_10h=0.0
    arfs_tb_10v=0.0
    arfs_tb_18h=0.0
    arfs_tb_18v=0.0
    arfs_tb_23h=0.0
    arfs_tb_23v=0.0
    arfs_tb_36h=0.0
    arfs_tb_36v=0.0
    arfs_tb_89h=0.0
    arfs_tb_89v=0.0
    arfs_land_water_frac=0.0
     !arfs_rfi_flag=0 ! uncomment for RFI_flag
     
     arfs_wt_tim=0.0
     arfs_wt_tb10v=0.0
     arfs_wt_tb10h=0.0
     arfs_wt_tb18v=0.0
     arfs_wt_tb18h=0.0
     arfs_wt_tb23v=0.0
     arfs_wt_tb23h=0.0
     arfs_wt_tb36v=0.0
     arfs_wt_tb36h=0.0
     arfs_wt_tb89v=0.0
     arfs_wt_tb89h=0.0
     arfs_wt_land_water_frac=0.0
     !arfs_wt_rfi_flag=0.0 ! uncomment for RFI flag
     
     arfs_quality_flag = 0
     snow_count = 0
     precip_count = 0
     ocean_count = 0
     total_count = 0
     excluded_snow_count = 0
     excluded_precip_count = 0

    ! Boundary check and debugging prints
    write(LDT_logunit,*) '[DEBUG] array dimensions:'
    write(LDT_logunit,*) '   nrows_l1rtb, ncols_l1rtb = ', nrows_l1rtb, ncols_l1rtb
    write(LDT_logunit,*) '   size(ref_lat), size(ref_lon) = ', size(ref_lat), size(ref_lon)


    DO ii = 1,ncols_l1rtb
        DO jj = 1,nrows_l1rtb
            ! Skip invalid coordinates
            if (lat_l1r(jj,ii) < -90.0 .or. lat_l1r(jj,ii) > 90.0 .or. &
                lon_l1r(jj,ii) < -180.0 .or. lon_l1r(jj,ii) > 180.0) then
                cycle
            endif
            
            ! FIND ARFS_GRID (r,c)
            c = MINLOC(ABS(lat_l1r(jj,ii)-ref_lat(:)),1) !Lat Direction
            r = MINLOC(ABS(lon_l1r(jj,ii)-ref_lon(:)),1) !Lon Direction
            
            ! Ensure r and c are valid indices
            if (r < 1 .or. r > size(ref_lon) .or. c < 1 .or. c > size(ref_lat)) then
                cycle  ! Skip this point
            endif
            
            rmin=r-5 ; IF (rmin < 1) rmin=1
            rmax=r+5 ; IF (rmax > size(ref_lon)) rmax=size(ref_lon)
            cmin=c-5 ; IF (cmin < 1) cmin=1
            cmax=c+5 ; IF (cmax > size(ref_lat)) cmax=size(ref_lat)
            
            ! Check snow and precip flags - with bounds protection
            IF (jj <= size(snow_flag,1) .and. ii <= size(snow_flag,2) .and. &
                jj <= size(precip_flag,1) .and. ii <= size(precip_flag,2)) THEN
                
                ! Determine if this footprint has snow or precip
                has_snow = (IBITS(snow_flag(jj,ii),qualitybit,1) == 1)
                has_precip = (IBITS(precip_flag(jj,ii),qualitybit,1) == 1)
                
                k=0
                DO rr = rmin,rmax !Lon direction
                    DO cc = cmin,cmax !Lat direction
                        lat1 = DBLE(lat_l1r(jj,ii)*d2r)
                        lon1 = DBLE(lon_l1r(jj,ii)*d2r)
                        lat2 = DBLE(ref_lat(cc)*d2r)
                        lon2 = DBLE(ref_lon(rr)*d2r)
                        
                        if(lat1.eq.lat2.and.lon1.eq.lon2) then
                            gcdist = 0.
                        else
                            gcdist = RE_KM * DACOS(DSIN(lat1) * DSIN(lat2) + DCOS(lat1) * DCOS(lat2) * DCOS(lon1-lon2))
                        endif
                        
                        IF (gcdist < search_radius) THEN !RESAMPLE ONLY WITHIN THE SEARCH RANGE
                            ! Always update counts for quality flag tracking (for ALL footprints)
                            total_count(rr,cc) = total_count(rr,cc) + 1
                            IF (IBITS(quality_flag(jj,ii), 0, 1) == 1) ocean_count(rr,cc) = ocean_count(rr,cc) + 1
                            IF (IBITS(quality_flag(jj,ii), 1, 1) == 1) precip_count(rr,cc) = precip_count(rr,cc) + 1  
                            IF (IBITS(quality_flag(jj,ii), 2, 1) == 1) snow_count(rr,cc) = snow_count(rr,cc) + 1
                            
                            ! Check if this footprint should be excluded
                            IF (has_snow .OR. has_precip) THEN
                                ! This footprint will be excluded from resampling
                                ! Just track that it was excluded (counts already updated above)
                                IF (has_snow) excluded_snow_count(rr,cc) = excluded_snow_count(rr,cc) + 1
                                IF (has_precip) excluded_precip_count(rr,cc) = excluded_precip_count(rr,cc) + 1
                            ELSE
                                ! No snow/precip - proceed with actual resampling
                                IF (gcdist < 0.0001D0) THEN !The TB is right on the grid center
                                    zerodistflag(rr,cc) = 1
                                    arfs_quality_flag(rr,cc) = quality_flag(jj,ii)
                                    
                                    IF ((ABS(tim(ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tim(rr,cc) = tim(ii)
                                        arfs_wt_tim(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_10h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_10h(rr,cc) = tb_10h(jj,ii)
                                        arfs_wt_tb10h(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_10v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_10v(rr,cc) = tb_10v(jj,ii)
                                        arfs_wt_tb10v(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_18h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_18h(rr,cc) = tb_18h(jj,ii)
                                        arfs_wt_tb18h(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_18v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_18v(rr,cc) = tb_18v(jj,ii)
                                        arfs_wt_tb18v(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_23h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_23h(rr,cc) = tb_23h(jj,ii)
                                        arfs_wt_tb23h(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_23v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_23v(rr,cc) = tb_23v(jj,ii)
                                        arfs_wt_tb23v(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_36h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_36h(rr,cc) = tb_36h(jj,ii)
                                        arfs_wt_tb36h(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_36v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_36v(rr,cc) = tb_36v(jj,ii)
                                        arfs_wt_tb36v(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_89h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_89h(rr,cc) = tb_89h(jj,ii)
                                        arfs_wt_tb89h(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(tb_89v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_tb_89v(rr,cc) = tb_89v(jj,ii)
                                        arfs_wt_tb89v(rr,cc) = 1.0
                                    ENDIF
                                    IF ((ABS(land_water_frac(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                        arfs_land_water_frac(rr,cc) = land_water_frac(jj,ii)
                                        arfs_wt_land_water_frac(rr,cc) = 1.0
                                    ENDIF
                                ELSE
                                    ! Weighted resampling (distance > 0)
                                    IF (zerodistflag(rr,cc).NE.1) THEN !Grid locations with an exact match are not updated
                                        k=k+1
                                        IF ((ABS(tim(ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tim(rr,cc) = arfs_tim(rr,cc) + tim(ii) * (1.0D0/gcdist)
                                            arfs_wt_tim(rr,cc) = arfs_wt_tim(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_10h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_10h(rr,cc) = arfs_tb_10h(rr,cc) + tb_10h(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb10h(rr,cc) = arfs_wt_tb10h(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_10v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_10v(rr,cc) = arfs_tb_10v(rr,cc) + tb_10v(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb10v(rr,cc) = arfs_wt_tb10v(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_18h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_18h(rr,cc) = arfs_tb_18h(rr,cc) + tb_18h(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb18h(rr,cc) = arfs_wt_tb18h(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_18v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_18v(rr,cc) = arfs_tb_18v(rr,cc) + tb_18v(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb18v(rr,cc) = arfs_wt_tb18v(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_23h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_23h(rr,cc) = arfs_tb_23h(rr,cc) + tb_23h(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb23h(rr,cc) = arfs_wt_tb23h(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_23v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_23v(rr,cc) = arfs_tb_23v(rr,cc) + tb_23v(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb23v(rr,cc) = arfs_wt_tb23v(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_36h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_36h(rr,cc) = arfs_tb_36h(rr,cc) + tb_36h(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb36h(rr,cc) = arfs_wt_tb36h(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_36v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_36v(rr,cc) = arfs_tb_36v(rr,cc) + tb_36v(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb36v(rr,cc) = arfs_wt_tb36v(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_89h(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_89h(rr,cc) = arfs_tb_89h(rr,cc) + tb_89h(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb89h(rr,cc) = arfs_wt_tb89h(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(tb_89v(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_tb_89v(rr,cc) = arfs_tb_89v(rr,cc) + tb_89v(jj,ii) * (1.0/gcdist)
                                            arfs_wt_tb89v(rr,cc) = arfs_wt_tb89v(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                        IF ((ABS(land_water_frac(jj,ii) - (-9999.0)).GT.1.0E-6)) THEN
                                            arfs_land_water_frac(rr,cc) = arfs_land_water_frac(rr,cc) + land_water_frac(jj,ii) * (1.0/gcdist)
                                            arfs_wt_land_water_frac(rr,cc) = arfs_wt_land_water_frac(rr,cc) + (1.0/gcdist)
                                        ENDIF
                                    ENDIF !zerodistflag check
                                ENDIF !(gcdist < 0.0001D0)
                            ENDIF !(has_snow .OR. has_precip)
                        ENDIF !(gcdist < search_radius)
                    END DO !cc = cmin,cmax
                END DO !rr = rmin,rmax
            ENDIF !(bounds check for snow_flag and precip_flag)
        END DO !jj=1,nrows_l1rtb
    END DO !ii=1,ncols_l1rtb

     ! TODO add a seperate for loop for the rfi_flag to loop trough lat89 and lon89 already defined in the variable defenition section but commented

    !APPLY WEIGHTING FUNCTION FOR RESAMPLING AND SET FILL VALUES
     WHERE(arfs_tim.NE.0.0 .AND. arfs_wt_tim.NE.0.0)
        arfs_tim = arfs_tim / arfs_wt_tim
     ELSEWHERE
        arfs_tim = -9999.0
     END WHERE
     
     WHERE(arfs_tb_10h.NE.0.0 .AND.arfs_wt_tb10h.NE.0.0)
        arfs_tb_10h= arfs_tb_10h / arfs_wt_tb10h
     ELSEWHERE
        arfs_tb_10h = -9999.0
     END WHERE
     
     WHERE(arfs_tb_10v.NE.0.0 .AND.arfs_wt_tb10v.NE.0.0)
        arfs_tb_10v= arfs_tb_10v / arfs_wt_tb10v
     ELSEWHERE
        arfs_tb_10v = -9999.0
     END WHERE
     
     WHERE(arfs_tb_18h.NE.0.0 .AND.arfs_wt_tb18h.NE.0.0)
        arfs_tb_18h= arfs_tb_18h / arfs_wt_tb18h
     ELSEWHERE
        arfs_tb_18h = -9999.0
     END WHERE
     
     WHERE(arfs_tb_18v.NE.0.0 .AND.arfs_wt_tb18v.NE.0.0)
        arfs_tb_18v= arfs_tb_18v / arfs_wt_tb18v
     ELSEWHERE
        arfs_tb_18v = -9999.0
     END WHERE
     
     WHERE(arfs_tb_23h.NE.0.0 .AND.arfs_wt_tb23h.NE.0.0)
        arfs_tb_23h= arfs_tb_23h / arfs_wt_tb23h
     ELSEWHERE
        arfs_tb_23h = -9999.0
     END WHERE
     
     WHERE(arfs_tb_23v.NE.0.0 .AND.arfs_wt_tb23v.NE.0.0)
        arfs_tb_23v= arfs_tb_23v / arfs_wt_tb23v
     ELSEWHERE
        arfs_tb_23v = -9999.0
     END WHERE
     
     WHERE(arfs_tb_36h.NE.0.0 .AND.arfs_wt_tb36h.NE.0.0)
        arfs_tb_36h= arfs_tb_36h / arfs_wt_tb36h
     ELSEWHERE
        arfs_tb_36h = -9999.0
     END WHERE
     
     WHERE(arfs_tb_36v.NE.0.0 .AND.arfs_wt_tb36v.NE.0.0)
        arfs_tb_36v= arfs_tb_36v / arfs_wt_tb36v
     ELSEWHERE
        arfs_tb_36v = -9999.0
     END WHERE
     
     WHERE(arfs_tb_89h.NE.0.0 .AND.arfs_wt_tb89h.NE.0.0)
        arfs_tb_89h= arfs_tb_89h / arfs_wt_tb89h
     ELSEWHERE
        arfs_tb_89h = -9999.0
     END WHERE
     
     WHERE(arfs_tb_89v.NE.0.0 .AND.arfs_wt_tb89v.NE.0.0)
        arfs_tb_89v= arfs_tb_89v / arfs_wt_tb89v
     ELSEWHERE
        arfs_tb_89v = -9999.0
     END WHERE

     WHERE(arfs_land_water_frac.NE.0.0 .AND.arfs_wt_land_water_frac.NE.0.0)
        arfs_land_water_frac = arfs_land_water_frac / arfs_wt_land_water_frac
     ELSEWHERE
        arfs_land_water_frac = -9999.0
     END WHERE

    ! Finalize quality flags using majority vote (only where no exact match occurred)
    DO i = 1, 2560
       DO j = 1, 1920
          IF (arfs_quality_flag(i,j) == 0) THEN
             ! Total footprints now already includes everything
             IF (total_count(i,j) > 0) THEN
                ! Ocean flag
                IF (ocean_count(i,j) > total_count(i,j)/2) THEN
                   arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 1)
                END IF
                
                ! Precipitation flag - note: excluded counts already in precip_count
                IF (precip_count(i,j) > total_count(i,j)/2) THEN
                   arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 2)
                END IF
                
                ! Snow flag - note: excluded counts already in snow_count
                IF (snow_count(i,j) > total_count(i,j)/2) THEN
                   arfs_quality_flag(i,j) = IOR(arfs_quality_flag(i,j), 4)
                END IF
             END IF
          END IF
       END DO
    END DO

     !WHERE(arfs_rfi_flag.NE.0.0 .AND.arfs_wt_rfi_flag.NE.0.0)
        !arfs_rfi_flag = arfs_rfi_flag / arfs_wt_rfi_flag
     !END WHERE
     ! Clean up allocated memory
     IF (allocated(zerodistflag)) DEALLOCATE(zerodistflag)
   END SUBROUTINE L1RTB2ARFS_INVDIS

   ! EMK...Only process subset of SMAP L1B fields for NRT operations.
   SUBROUTINE L1BTB2ARFS_INVDIS_SUBSET(tim, tbvl1b_cor, &
        lat_l1b, lon_l1b, tbv_qual_flag, tbh_qual_flag, &
        sc_nadir_angle, antenna_scan_angle, nrows_l1rtb, ncols_l1rtb, &
        ref_lat, ref_lon, arfs_tim, arfs_tbv_cor)

     use LDT_logMod, only: LDT_logunit

     INTEGER(4) :: ii, jj, k, r, c, rr, rmin, rmax, cc, cmin, cmax, nrows_l1rtb, ncols_l1rtb
     INTEGER(4), PARAMETER :: qualitybit = 0
     REAL(8), PARAMETER :: RE_KM = 6371.228, search_radius = 20.0, PI = 3.141592653589793238, d2r = PI/180.0
     REAL(8)  :: gcdist, lat1, lon1, lat2, lon2
     REAL*8,DIMENSION(nrows_l1rtb) :: tim
     REAL*4,DIMENSION(nrows_l1rtb,ncols_l1rtb) :: tbvl1b_cor
     REAL*4,DIMENSION(nrows_l1rtb,ncols_l1rtb) :: lat_l1b, lon_l1b, antenna_scan_angle
     REAL*4,DIMENSION(ncols_l1rtb) :: sc_nadir_angle
     INTEGER*4,DIMENSION(nrows_l1rtb,ncols_l1rtb) :: tbv_qual_flag, tbh_qual_flag
     INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: zerodistflag
     REAL*8,DIMENSION(:),ALLOCATABLE :: ref_lat, ref_lon
     REAL*8,DIMENSION(2560,1920) :: arfs_tim
     REAL*4,DIMENSION(2560,1920) :: arfs_tbv_cor
     REAL*8,DIMENSION(2560,1920) :: arfs_wt_tim
     REAL*4,DIMENSION(2560,1920) :: arfs_wt_cor_tbv
     INTEGER*4,DIMENSION(2560,1920) :: arfs_samplenumv

     ALLOCATE(zerodistflag(size(ref_lon),size(ref_lat)))
     zerodistflag = 0 ! EMK Added initialization
     
     !INITIAL THE OUTPUT VARIABLES
     arfs_tim=0.0
     arfs_tbv_cor=0.0
     arfs_wt_tim=0.0
     arfs_wt_cor_tbv=0.0
     arfs_samplenumv=0.0

     DO ii = 1,ncols_l1rtb
        IF (ABS (sc_nadir_angle(ii)) <= 2.0) THEN
           DO jj = 1,nrows_l1rtb

              
              IF (ABS (antenna_scan_angle(jj,ii)).LE.360.00) THEN

                 lat1 = DBLE (lat_l1b(jj,ii)*d2r)
                 lon1 = DBLE (lon_l1b(jj,ii)*d2r)
                 ! FIND ARFS_GRID (r,c)
                 c = MINLOC(ABS(lat_l1b(jj,ii)-ref_lat(:)),1) !Lat Direction
                 r = MINLOC(ABS(lon_l1b(jj,ii)-ref_lon(:)),1) !Lon Direction
                 rmin=r-5 ; IF (rmin < 1) rmin=1
                 rmax=r+5 ; IF (rmax > size(ref_lon)) rmax=size(ref_lon)
                 cmin=c-5 ; IF (cmin < 1) cmin=1
                 cmax=c+5 ; IF (cmax > size(ref_lat)) cmax=size(ref_lat)
!                 IF (IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(jj,jj),qualitybit,1) == 0) THEN !RESAMPLE ONLY WHEN BOTH V and H MEET QUALITY
                 IF (IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(jj,ii),qualitybit,1) == 0) THEN !RESAMPLE ONLY WHEN BOTH V and H MEET QUALITY

                    k=0
                    !                         DO rr = rmin,rmax !Lon direction
                    !                            DO cc =cmin,cmax !Lat direction
                    DO cc =cmin,cmax !Lat direction
                       lat2 = DBLE (ref_lat(cc)*d2r)
                       DO rr = rmin,rmax !Lon direction
                          lon2 = DBLE (ref_lon(rr)*d2r)

                          if(lat1.eq.lat2.and.lon1.eq.lon2) then
                             gcdist = 0.
                          else
                             gcdist = RE_KM * DACOS ( DSIN (lat1) * DSIN (lat2) + DCOS (lat1) * DCOS (lat2) * DCOS (lon1-lon2) )
                          endif

                          IF (gcdist < search_radius) THEN !RESAMPLE ONLY WITHIN THE SEARCH RANGE
                             IF (gcdist < 0.0001D0) THEN !The TB is right on the grid center
                                zerodistflag (rr,cc) = 1
                                !IF ((ABS (tim(ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                IF ( .not. tim(ii) < 0) THEN !DO IF NOT FILLVALUE(-9999)

                                   arfs_tim(rr,cc) = tim(ii) ; arfs_wt_tim(rr,cc) = 1.0
                                END IF
                                !IF ((ABS (tbvl1b_cor(jj,ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                IF ( .not. tbvl1b_cor(jj,ii) < 0 ) THEN !DO IF NOT FILLVALUE(-9999)

                                   arfs_tbv_cor(rr,cc) = tbvl1b_cor(jj,ii) ; arfs_wt_cor_tbv(rr,cc) = 1.0
                                   arfs_samplenumv(rr,cc)=1 !Sample number only calculate for correct tb
                                   k=k+1;
                                END IF
                             ELSE
                                IF (zerodistflag (rr,cc).EQ.0) THEN

                                   !IF ((ABS (tim(ii) - (-9999.0)).GT.1.0D-7)) THEN !DO IF NOT FILLVALUE(-9999)
                                   IF ( .not. tim(ii) < 0 ) THEN !DO IF NOT FILLVALUE(-9999)

                                      !arfs_tim(rr,cc) = arfs_tim(rr,cc) + tim(ii) / SNGL (gcdist*gcdist)
                                      !arfs_wt_tim(rr,cc) = arfs_wt_tim(rr,cc) + 1.0 / SNGL (gcdist*gcdist)
                                      arfs_tim(rr,cc) = arfs_tim(rr,cc) + tim(ii) / (gcdist*gcdist)
                                      arfs_wt_tim(rr,cc) = arfs_wt_tim(rr,cc) + 1.0 /  (gcdist*gcdist)

                                   END IF
                                   !IF ((ABS (tbvl1b_cor(jj,ii) - (-9999.0))).GT.1.0D-7) THEN !DO IF NOT FILLVALUE(-9999)
                                  IF ( .not. tbvl1b_cor(jj,ii) < 0 ) THEN !DO IF NOT FILLVALUE(-9999)
 
                                      arfs_tbv_cor(rr,cc) = arfs_tbv_cor(rr,cc) + tbvl1b_cor(jj,ii) / SNGL (gcdist*gcdist)
                                      arfs_wt_cor_tbv(rr,cc) = arfs_wt_cor_tbv(rr,cc) + 1.0 /  SNGL (gcdist*gcdist)
                                      arfs_samplenumv(rr,cc)=arfs_samplenumv(rr,cc)+1.0 !Sample number only calculate for correct tb
                                      k=k+1;
                                   END IF

                                END IF !(zerodistflag (rr,cc) = 0)
                             END IF !(gcdist < 0.0001D0)!
                          END IF !(gcdist < search_radius)
                          !                            END DO !cc =cmin,cmax
                          !                         END DO !rr = rmin,rmax
                       END DO !rr = rmin,rmax
                    END DO !cc =cmin,cmax

                 END IF !(IBITS (tbv_qual_flag(jj,ii),qualitybit,1) == 0 .AND. IBITS (tbh_qual_flag(ii,jj),qualitybit,1) == 0)
              END IF !(ABS (antenna_scan_angle(ii,jj)) <= 360.00)
           END DO !jj=1,2
        END IF !(ABS (sc_nadir_angle(ii)) <= 2.0)
     END DO !ii=1,2

     !APPLY WEIGHTING FUNCTION FOR RESAMPLING
     WHERE(arfs_tim.NE.0.0.AND.arfs_wt_tim.NE.0.0)
        arfs_tim = arfs_tim / arfs_wt_tim
     END WHERE
     WHERE(arfs_tbv_cor.NE.0.0.AND.arfs_wt_cor_tbv.NE.0.0)
        arfs_tbv_cor = arfs_tbv_cor / arfs_wt_cor_tbv
     END WHERE

     deallocate(zerodistflag) ! EMK cleanup memory
   END SUBROUTINE L1BTB2ARFS_INVDIS_SUBSET

 END MODULE invdist_l1r2arfs
