!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! MODULE: TOOLSUBS_AMSR3
!
! DESCRIPTION: Reads AMSR3 L1R (netCDF-4) brightness temperatures and their
!   per-channel quality onto the native 243-point P89o grid. TB and quality
!   are returned as 3D arrays (nfov, nscan, nchan) driven by a channel table,
!   so adding/removing a channel is a one-row edit. Snow/precip are computed
!   per footprint (the only masks applied before resampling); the raw
!   per-channel quality bytes are carried through untouched.
!-------------------------------------------------------------------------

#include "LDT_misc.h"

MODULE TOOLSUBS_AMSR3
    USE LDT_logMod, only: LDT_logunit, LDT_endrun
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    USE netcdf
#endif
    IMPLICIT NONE

    INTEGER, PARAMETER :: NCH = 14   ! channels read (edit the table below to change)

CONTAINS

    SUBROUTINE get_amsr3_l1r(filename, tb, tbq, ch_name, ch_freq, ch_pol, &
         lat, lon, tb_time, land_frac, snow, precip, scan_qual, nfov, nscan, nchan, ierr)

      character(*), intent(in) :: filename
      real*4,       allocatable, intent(out) :: tb(:,:,:)      ! (nfov,nscan,nchan) K
      integer*2,    allocatable, intent(out) :: tbq(:,:,:)     ! (nfov,nscan,nchan) raw quality byte (0-255)
      integer*2,    allocatable, intent(out) :: scan_qual(:)   ! (nscan) ScanDataQuality byte; nonzero = bad scan
      character(len=20), intent(out) :: ch_name(NCH)
      real*4,            intent(out) :: ch_freq(NCH)
      character(len=1),  intent(out) :: ch_pol(NCH)
      real*4,    allocatable, intent(out) :: lat(:,:), lon(:,:)
      real*8,    allocatable, intent(out) :: tb_time(:)        ! (nscan) TAI93 seconds
      integer*4, allocatable, intent(out) :: land_frac(:,:), snow(:,:), precip(:,:)
      integer,   intent(out) :: nfov, nscan, nchan, ierr

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
      integer :: ncid, did, k, i, j
      integer :: i18v, i18h, i23v, i36v, i89v
      logical :: exists
      real    :: sil, tt18

      ierr = 1 ; nfov = 0 ; nscan = 0 ; nchan = NCH

      ! ---- channel table: name | frequency (GHz) | polarization -------------
      ! native FOV per band: 6.9/7.3 -> FOV06, 10.65 -> FOV10, 18.7/23.8 -> FOV23,
      ! 36.5/89 -> FOV36. Add a row (and bump NCH) to include e.g. 10.25 (Ch10u).
      ch_name(1) ='Tb_FOV06Ch06V_P89o' ; ch_freq(1) = 6.9  ; ch_pol(1) ='V'
      ch_name(2) ='Tb_FOV06Ch06H_P89o' ; ch_freq(2) = 6.9  ; ch_pol(2) ='H'
      ch_name(3) ='Tb_FOV06Ch07V_P89o' ; ch_freq(3) = 7.3  ; ch_pol(3) ='V'
      ch_name(4) ='Tb_FOV06Ch07H_P89o' ; ch_freq(4) = 7.3  ; ch_pol(4) ='H'
      ch_name(5) ='Tb_FOV10Ch10V_P89o' ; ch_freq(5) =10.65 ; ch_pol(5) ='V'
      ch_name(6) ='Tb_FOV10Ch10H_P89o' ; ch_freq(6) =10.65 ; ch_pol(6) ='H'
      ch_name(7) ='Tb_FOV23Ch18V_P89o' ; ch_freq(7) =18.7  ; ch_pol(7) ='V'
      ch_name(8) ='Tb_FOV23Ch18H_P89o' ; ch_freq(8) =18.7  ; ch_pol(8) ='H'
      ch_name(9) ='Tb_FOV23Ch23V_P89o' ; ch_freq(9) =23.8  ; ch_pol(9) ='V'
      ch_name(10)='Tb_FOV23Ch23H_P89o' ; ch_freq(10)=23.8  ; ch_pol(10)='H'
      ch_name(11)='Tb_FOV36Ch36V_P89o' ; ch_freq(11)=36.5  ; ch_pol(11)='V'
      ch_name(12)='Tb_FOV36Ch36H_P89o' ; ch_freq(12)=36.5  ; ch_pol(12)='H'
      ch_name(13)='Tb_FOV36Ch89V_P89o' ; ch_freq(13)=89.0  ; ch_pol(13)='V'
      ch_name(14)='Tb_FOV36Ch89H_P89o' ; ch_freq(14)=89.0  ; ch_pol(14)='H'

      inquire(file=trim(filename), exist=exists)
      if (.not. exists) then
         write(LDT_logunit,*)'[ERR] Cannot find file ', trim(filename) ; return
      end if

      ! ---- open + learn swath size (nfov = OBS_89o pixels, nscan = time) -----
      if (.not. ok(nf90_open(trim(filename), nf90_nowrite, ncid), 'open')) return
      if (.not. ok(nf90_inq_dimid(ncid,'OBS_89o',did),'OBS_89o')) return
      if (.not. ok(nf90_inquire_dimension(ncid,did,len=nfov),'OBS_89o len')) return
      if (.not. ok(nf90_inq_dimid(ncid,'time',did),'time')) return
      if (.not. ok(nf90_inquire_dimension(ncid,did,len=nscan),'time len')) return
      write(LDT_logunit,*)'[INFO] AMSR3 grid nfov=',nfov,' nscan=',nscan,' nchan=',NCH

      allocate(tb(nfov,nscan,NCH), tbq(nfov,nscan,NCH), lat(nfov,nscan), lon(nfov,nscan), &
               tb_time(nscan), land_frac(nfov,nscan), snow(nfov,nscan), precip(nfov,nscan), &
               scan_qual(nscan))
      snow = 0 ; precip = 0

      ! ---- read every channel + its quality in one loop ---------------------
      do k = 1, NCH
         if (.not. rd_tb(ch_name(k), tb(:,:,k))) return
         if (.not. rd_q (trim(ch_name(k))//'_Quality', tbq(:,:,k))) return
      end do

      ! ---- geolocation, time, land fraction ---------------------------------
      if (.not. rd_real('Latitude_P89o',  lat)) return
      if (.not. rd_real('Longitude_P89o', lon)) return
      if (.not. rd_time(tb_time)) return
      if (.not. rd_int('LandAreaPercent_FOV06_P89o', land_frac)) return
      if (.not. rd_q1d('ScanDataQuality', scan_qual)) return   ! bits 3-7: missing/nav/attitude/HTS/antenna

      call ok_warn(nf90_close(ncid), 'close')

      ! ---- snow / precip per footprint (only masks applied pre-resampling) ---
      i18v=chan_idx(18.7,'V') ; i18h=chan_idx(18.7,'H') ; i23v=chan_idx(23.8,'V')
      i36v=chan_idx(36.5,'V') ; i89v=chan_idx(89.0,'V')
      if (min(i18v,i18h,i23v,i36v,i89v) > 0) then
         do j = 1, nscan
            do i = 1, nfov
               if (land_frac(i,j) >= 50 .and. tb(i,j,i18v) > 0 .and. tb(i,j,i18h) > 0 .and. &
                   tb(i,j,i23v) > 0 .and. tb(i,j,i36v) > 0 .and. tb(i,j,i89v) > 0) then
                  sil  = 451.88 - 0.44*tb(i,j,i18v) - 1.775*tb(i,j,i23v) + &
                         0.00574*tb(i,j,i23v)**2 - tb(i,j,i89v)
                  tt18 = tb(i,j,i18v) - tb(i,j,i18h)
                  if (sil > 10) then
                     if (tb(i,j,i23v) <= 264.0 .and. tb(i,j,i23v) <= (175.0 + 0.49*tb(i,j,i89v))) then
                        snow(i,j) = 1
                        if (tt18 >= 18 .and. (tb(i,j,i18v)-tb(i,j,i36v)) <= 10 .and. &
                            (tb(i,j,i36v)-tb(i,j,i89v)) <= 10) snow(i,j) = 0
                        if (tt18 >= 8  .and. (tb(i,j,i18v)-tb(i,j,i36v)) <= 2  .and. &
                            (tb(i,j,i23v)-tb(i,j,i89v)) <= 6)  snow(i,j) = 1
                     else
                        precip(i,j) = 1
                        if (tt18 > 20) precip(i,j) = 0
                        if (tb(i,j,i89v) > 253 .and. tt18 > 7) precip(i,j) = 0
                     endif
                  endif
               endif
            end do
         end do
      endif

      ierr = 0

    CONTAINS

      logical function ok(status, tag)
        integer, intent(in) :: status ; character(*), intent(in) :: tag
        ok = (status == nf90_noerr)
        if (.not. ok) write(LDT_logunit,*)'[ERR] netCDF ', trim(tag), ': ', trim(nf90_strerror(status))
      end function ok

      subroutine ok_warn(status, tag)
        integer, intent(in) :: status ; character(*), intent(in) :: tag
        if (status /= nf90_noerr) write(LDT_logunit,*)'[WARN] netCDF ', trim(tag)
      end subroutine ok_warn

      ! index of the channel with this frequency and polarization (0 if absent)
      integer function chan_idx(f, p)
        real, intent(in) :: f ; character(len=1), intent(in) :: p
        integer :: kk
        chan_idx = 0
        do kk = 1, NCH
           if (abs(ch_freq(kk)-f) < 0.5 .and. ch_pol(kk) == p) then
              chan_idx = kk ; return
           endif
        end do
      end function chan_idx

      ! uint16 * 0.01; 65534 (missing) and 65535 (parity) -> -9999. Read into
      ! int32 so unsigned values above 32767 do not wrap negative.
      logical function rd_tb(name, out2d)
        character(*), intent(in) :: name ; real*4, intent(out) :: out2d(nfov,nscan)
        integer :: vid ; integer*4, allocatable :: raw(:,:)
        rd_tb = .false.
        if (.not. ok(nf90_inq_varid(ncid, name, vid), name)) return
        allocate(raw(nfov,nscan))
        if (.not. ok(nf90_get_var(ncid, vid, raw), name)) return
        where (raw >= 65534) ; out2d = -9999.0 ; elsewhere ; out2d = real(raw)*0.01 ; end where
        deallocate(raw) ; rd_tb = .true.
      end function rd_tb

      ! quality byte: read uint8 into int16 (values reach 128/255, past int8 range)
      logical function rd_q(name, out2d)
        character(*), intent(in) :: name ; integer*2, intent(out) :: out2d(nfov,nscan)
        integer :: vid ; rd_q = .false.
        if (.not. ok(nf90_inq_varid(ncid, name, vid), name)) return
        if (.not. ok(nf90_get_var(ncid, vid, out2d), name)) return
        rd_q = .true.
      end function rd_q

      logical function rd_q1d(name, out1d)      ! per-scan uint8 quality -> int16
        character(*), intent(in) :: name ; integer*2, intent(out) :: out1d(nscan)
        integer :: vid ; rd_q1d = .false.
        if (.not. ok(nf90_inq_varid(ncid, name, vid), name)) return
        if (.not. ok(nf90_get_var(ncid, vid, out1d), name)) return
        rd_q1d = .true.
      end function rd_q1d

      logical function rd_int(name, out2d)      ! uint8 land percent -> int32
        character(*), intent(in) :: name ; integer*4, intent(out) :: out2d(nfov,nscan)
        integer :: vid ; rd_int = .false.
        if (.not. ok(nf90_inq_varid(ncid, name, vid), name)) return
        if (.not. ok(nf90_get_var(ncid, vid, out2d), name)) return
        rd_int = .true.
      end function rd_int

      logical function rd_real(name, out2d)     ! float32 lat/lon, fill -9999
        character(*), intent(in) :: name ; real*4, intent(out) :: out2d(nfov,nscan)
        integer :: vid ; rd_real = .false.
        if (.not. ok(nf90_inq_varid(ncid, name, vid), name)) return
        if (.not. ok(nf90_get_var(ncid, vid, out2d), name)) return
        where (out2d < -9000.0) out2d = -9999.0
        rd_real = .true.
      end function rd_real

      logical function rd_time(out1d)           ! float64 seconds since 1993-01-01 (TAI)
        real*8, intent(out) :: out1d(nscan) ; integer :: vid ; rd_time = .false.
        if (.not. ok(nf90_inq_varid(ncid,'ScanTimeTAI93',vid),'ScanTimeTAI93')) return
        if (.not. ok(nf90_get_var(ncid, vid, out1d),'ScanTimeTAI93')) return
        rd_time = .true.
      end function rd_time

#else
      write(LDT_logunit,*)'[ERR] get_amsr3_l1r called without netCDF support!'
      call LDT_endrun()
      ierr = 1
#endif

    END SUBROUTINE get_amsr3_l1r
END MODULE TOOLSUBS_AMSR3
