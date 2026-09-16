!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: invdist_l1r2arfs_AMSR3
!
! DESCRIPTION: Inverse-distance resampling of AMSR3 L1R swath TB (all
!   channels at once, 3D) onto the Air Force (ARFS) grid.
!
!   EXCLUDED before accumulating (validity):
!     footprint : bad scan (ScanDataQuality /= 0), bad lat/lon,
!                 geometry abnormal (quality bit 2)
!     per chan  : TB fill, TB abnormal (bit 3), resampling abnormal (bits 6-5 == 11)
!   CONFIGURABLE (filter_snow_precip): snow/precip footprints are dropped
!     from the TB average when .true.; kept when .false. (snow products).
!     Their weighted share per cell is carried either way.
!   PER-CHANNEL radius_km(k): each channel collects footprints within its
!     own radius, so fine channels (36, 89 GHz) are not over-smoothed by a
!     radius chosen for the 6.9 GHz footprint.
!   CARRIED through, weighted by the same 1/dist used for TB:
!     arfs_frfi : share with RFI occurred (bits 1-0 == 10)
!     arfs_fcau : share with resampling caution (bits 6-5 == 10)
!     arfs_qor  : bitwise OR of every contributing quality byte
!     arfs_fsnow, arfs_fprec : share of footprints flagged snow / precip
!-------------------------------------------------------------------------

MODULE invdist_l1r2arfs_AMSR3
  USE LDT_logMod, only: LDT_logunit
  IMPLICIT NONE
CONTAINS

  SUBROUTINE L1RTB2ARFS_INVDIS_AMSR3(tb, tbq, snow, precip, scan_qual, use_scan, land_frac, tim, &
       lat, lon, nfov, nscan, nchan, ref_lat, ref_lon, radius_km, filter_snow_precip, &
       arfs_tb, arfs_frfi, arfs_fcau, arfs_qor, arfs_fsnow, arfs_fprec, &
       arfs_tim, arfs_land, arfs_nsamp)

    integer,   intent(in) :: nfov, nscan, nchan
    real*4,    intent(in) :: tb(nfov,nscan,nchan), lat(nfov,nscan), lon(nfov,nscan)
    integer*2, intent(in) :: tbq(nfov,nscan,nchan), scan_qual(nscan)
    logical,   intent(in) :: use_scan(nscan)                ! .false. = scan belongs to the other pass
    integer*4, intent(in) :: snow(nfov,nscan), precip(nfov,nscan), land_frac(nfov,nscan)
    real*8,    intent(in) :: tim(nscan), ref_lat(:), ref_lon(:)
    real*4,    intent(in) :: radius_km(nchan)          ! search radius per channel
    logical,   intent(in) :: filter_snow_precip        ! .true. = drop snow/precip footprints from TB
    real*4,    allocatable, intent(out) :: arfs_tb(:,:,:), arfs_frfi(:,:,:), arfs_fcau(:,:,:)
    real*4,    allocatable, intent(out) :: arfs_fsnow(:,:), arfs_fprec(:,:), arfs_land(:,:)
    integer*2, allocatable, intent(out) :: arfs_qor(:,:,:)
    real*8,    allocatable, intent(out) :: arfs_tim(:,:)
    integer*4, allocatable, intent(out) :: arfs_nsamp(:,:)

    real*8, parameter :: RE_KM = 6371.228d0, d2r = 3.141592653589793d0/180.d0
    integer   :: nlon, nlat, ii, jj, k, r, c, rr, cc, nlatw, nlonw
    integer   :: nskip_scan, nskip_geo, nskip_sp
    integer*2 :: q
    logical   :: sp
    real*8    :: lat1, lon1, lat2, lon2, gcdist, w, rmax, dlat_km, dlon_km
    real*4, allocatable :: wt(:,:,:), wt_fp(:,:), wt_tim(:,:), wt_land(:,:)

    nlon = size(ref_lon) ; nlat = size(ref_lat)
    rmax    = maxval(radius_km)
    dlat_km = abs(ref_lat(2)-ref_lat(1)) * 111.2d0          ! ARFS cell size in km, from the grid itself
    dlon_km = abs(ref_lon(2)-ref_lon(1)) * 111.2d0          ! at the equator; shrinks by cos(lat)
    nlatw   = ceiling(rmax/dlat_km) + 1                     ! window half-width in cells (lat)

    allocate(arfs_tb(nlon,nlat,nchan), arfs_frfi(nlon,nlat,nchan), arfs_fcau(nlon,nlat,nchan), &
             arfs_qor(nlon,nlat,nchan), wt(nlon,nlat,nchan), &
             arfs_fsnow(nlon,nlat), arfs_fprec(nlon,nlat), wt_fp(nlon,nlat), &
             arfs_tim(nlon,nlat), arfs_land(nlon,nlat), arfs_nsamp(nlon,nlat), wt_tim(nlon,nlat), wt_land(nlon,nlat))
    arfs_tb = 0 ; arfs_frfi = 0 ; arfs_fcau = 0 ; arfs_qor = 0 ; wt = 0     ! fractions accumulate in place
    arfs_fsnow = 0 ; arfs_fprec = 0 ; wt_fp = 0
    arfs_tim = 0 ; arfs_land = 0 ; arfs_nsamp = 0 ; wt_tim = 0 ; wt_land = 0
    nskip_scan = 0 ; nskip_geo = 0 ; nskip_sp = 0

    do ii = 1, nscan
       if (.not. use_scan(ii)) cycle                  ! other pass (asc/desc); not an error, not counted
       if (scan_qual(ii) /= 0) then                    ! bits 0-2 fixed 0, so any nonzero = bit 3-7 error or fill 255
          nskip_scan = nskip_scan + 1 ; cycle
       endif
       do jj = 1, nfov
          ! ---- footprint-level validity exclusions ---------------------------
          lat1 = lat(jj,ii) ; lon1 = lon(jj,ii)
          if (lon1 > 180.d0) lon1 = lon1 - 360.d0
          if (lat1 < -90.d0 .or. lat1 > 90.d0 .or. lon1 < -180.d0 .or. lon1 > 180.d0) cycle
          if (any(ibits(tbq(jj,ii,:), 2, 1) == 1)) then                     ! geometry abnormal
             nskip_geo = nskip_geo + 1 ; cycle
          endif
          sp = (snow(jj,ii) == 1 .or. precip(jj,ii) == 1)                  ! counted below, dropped only if configured
          if (filter_snow_precip .and. sp) nskip_sp = nskip_sp + 1

          ! ---- nearest cell, window sized for the largest radius -------------
          c = minloc(abs(lat1 - ref_lat), 1)
          r = minloc(abs(lon1 - ref_lon), 1)
          nlonw = ceiling(rmax / max(dlon_km*cos(d2r*lat1), 0.5d0)) + 1   ! lon cells narrow toward the poles

          do rr = max(1, r-nlonw), min(nlon, r+nlonw)
             do cc = max(1, c-nlatw), min(nlat, c+nlatw)
                lat2 = ref_lat(cc) ; lon2 = ref_lon(rr)
                gcdist = 2.d0*RE_KM*asin(min(1.d0, sqrt(sin(d2r*(lat1-lat2)/2)**2 + &
                         cos(d2r*lat1)*cos(d2r*lat2)*sin(d2r*(lon1-lon2)/2)**2)))
                if (gcdist >= rmax) cycle
                w = 1.d0 / max(gcdist, 1.d-3)

                ! ---- footprint-level layers (all footprints in the window) -------
                arfs_nsamp(rr,cc) = arfs_nsamp(rr,cc) + 1
                wt_fp(rr,cc)      = wt_fp(rr,cc) + w
                if (snow(jj,ii)   == 1) arfs_fsnow(rr,cc) = arfs_fsnow(rr,cc) + w
                if (precip(jj,ii) == 1) arfs_fprec(rr,cc) = arfs_fprec(rr,cc) + w
                if (tim(ii) > -9000.d0) then
                   arfs_tim(rr,cc) = arfs_tim(rr,cc) + tim(ii)*w ; wt_tim(rr,cc) = wt_tim(rr,cc) + w
                endif
                if (land_frac(jj,ii) >= 0 .and. land_frac(jj,ii) <= 100) then
                   arfs_land(rr,cc) = arfs_land(rr,cc) + land_frac(jj,ii)*w ; wt_land(rr,cc) = wt_land(rr,cc) + w
                endif
                if (filter_snow_precip .and. sp) cycle                   ! keep the layers above, skip the TB

                ! ---- per channel: own radius, exclude on validity, carry the rest ---
                do k = 1, nchan
                   if (gcdist >= radius_km(k)) cycle                     ! outside this channel's radius
                   q = tbq(jj,ii,k)
                   if (tb(jj,ii,k) < -9000.0) cycle                      ! fill
                   if (ibits(q,3,1) == 1)     cycle                      ! TB abnormal
                   if (ibits(q,5,2) == 3)     cycle                      ! resampling abnormal (11)
                   arfs_tb(rr,cc,k) = arfs_tb(rr,cc,k) + tb(jj,ii,k)*w
                   wt(rr,cc,k)      = wt(rr,cc,k) + w
                   if (ibits(q,0,2) == 2) arfs_frfi(rr,cc,k) = arfs_frfi(rr,cc,k) + w   ! RFI occurred (10)
                   if (ibits(q,5,2) == 2) arfs_fcau(rr,cc,k) = arfs_fcau(rr,cc,k) + w   ! resampling caution (10)
                   arfs_qor(rr,cc,k) = ior(arfs_qor(rr,cc,k), q)
                end do
             end do
          end do
       end do
    end do

    ! ---- normalize by the same weights ----------------------------------------
    where (wt > 0)
       arfs_tb   = arfs_tb   / wt
       arfs_frfi = arfs_frfi / wt
       arfs_fcau = arfs_fcau / wt
    elsewhere
       arfs_tb = -9999.0 ; arfs_frfi = -9999.0 ; arfs_fcau = -9999.0 ; arfs_qor = -1
    end where
    where (wt_fp   > 0) ; arfs_fsnow = arfs_fsnow / wt_fp ; arfs_fprec = arfs_fprec / wt_fp
    elsewhere           ; arfs_fsnow = -9999.0 ; arfs_fprec = -9999.0 ; end where
    where (wt_tim  > 0) ; arfs_tim  = arfs_tim  / wt_tim  ; elsewhere ; arfs_tim  = -9999.d0 ; end where
    where (wt_land > 0) ; arfs_land = arfs_land / wt_land ; elsewhere ; arfs_land = -9999.0  ; end where

    write(LDT_logunit,*)'[INFO] AMSR3 resample: radius km min/max=', minval(radius_km), maxval(radius_km), &
         ' filter_snow_precip=', filter_snow_precip
    write(LDT_logunit,*)'[INFO] skipped scans=', nskip_scan, ' geom-bad=', nskip_geo, &
         ' snow/precip dropped=', nskip_sp, ' | cells filled=', count(arfs_nsamp > 0)
    deallocate(wt, wt_fp, wt_tim, wt_land)

  END SUBROUTINE L1RTB2ARFS_INVDIS_AMSR3
END MODULE invdist_l1r2arfs_AMSR3
