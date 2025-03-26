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
!  25 Dec 2024: Ehsan Jalilvand; AMSR2 implementation
!
! DESCRIPTION: RESAMPLE AMSR_L1R TB TO AIR FORCE GRID
! INPUT : AMSR - L1R Brightness Temperature
! OUTPUT: AMSRTB_ARFSGRIDE_ddmmyyy.dat
! NOTES : Inverse Distance Squared with 0.4 deg serching window
!-------------------------------------------------------------------------

subroutine AMSR_L1R_SUBSET(AMSRFILE,L1R_dir,Orbit,ARFS_TIME,rc)

  USE VARIABLES
  USE DATADOMAIN
  USE FUNCTIONS
  USE TOOLSUBS
  USE invdist_l1r2arfs
  USE LDT_logMod
  USE LDT_smap_e_oplMod

  IMPLICIT NONE

  INTEGER :: i, j, nrow, mcol
  CHARACTER (len=100) :: AMSRFILE
  character (len=100) :: L1R_dir
  character (len=20)  :: variable_name(10)
  character (len=100) :: resample_filename(10)
  character (len=1)   :: Orbit ! check where these files are being used
  integer             :: var_i
  integer             :: L1R_dir_len,L1R_fname_len
  integer :: ierr
  integer :: rc

  REAL*8,DIMENSION(:,:),ALLOCATABLE :: TIME_L1R
  REAL*4,DIMENSION(:,:),ALLOCATABLE :: TB_10H, TB_10V, TB_18H, TB_18V, TB_23H, TB_23V, TB_36H, TB_36V, TB_89H, TB_89V
  REAL*4,DIMENSION(:,:),ALLOCATABLE :: LAT_L1R, LON_L1R
! REAL*4,DIMENSION(:),ALLOCATABLE :: ANTSCN_L1B, SCNANG_L1B !Which AMSR layer should be used instead of SMAP anthena scan angle and
! scan angle
  INTEGER*4,DIMENSION(:,:),ALLOCATABLE :: FLAG_L1R
  
  REAL*8,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
  REAL*8,DIMENSION(2560,1920) :: ARFS_TIME
  REAL*4,DIMENSION(2560,1920) :: ARFS_TB_10V, ARFS_TB_10H, ARFS_TB_18H, ARFS_TB_18V,ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, ARFS_TB_89H, ARFS_TB_89V ! This is instead of ARFS_COR_TBV
  INTEGER*4,DIMENSION(2560,1920),ALLOCATABLE :: ARFS_FLAG_L1R

  REAL :: T1, T2

  rc = 0

  CALL ARFS_GEO
  ALLOCATE(ARFS_LAT(arfs_nrow_lat),ARFS_LON(arfs_mcol_lon)) 
  ARFS_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_space)
  ARFS_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_space)


  CALL GetAMSR_L1R (AMSRFILE,TIME_L1R &
          TB_10H, TB_10V, TB_18H, TB_18V, &
          TB_23H, TB_23V, TB_36H, TB_36V, &
          TB_89H, TB_89V, LAT_L1R, LON_L1R, &
          FLAG_L1R, nrow, mcol, ierr)
 
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

  CALL L1RTB2ARFS_INVDIS_SUBSET(TIME_L1R, &
          TB_10H, TB_10V, TB_18H, TB_18V, &
          TB_23H, TB_23V, TB_36H, TB_36V, &
          TB_89H, TB_89V, &
          LAT_L1R, LON_L1R, FLAG_L1R, &
          nrow, mcol, &
          ARFS_LAT, ARFS_LON, ARFS_TIME, &
          ARFS_TB_10H, ARFS_TB_10V, ARFS_TB_18H, ARFS_TB_18V, &
          ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, &
          ARFS_TB_89H, ARFS_TB_89V, ARFS_FLAG_L1R)
  
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
  AMSReOPL%ARFS_FLAG_L1R = ARFS_FLAG_L1R

  L1R_dir_len = len_trim(L1R_dir)
  L1R_fname_len = len_trim(AMSRFILE)

  if(AMSReOPL%L1Rtype.eq.1) then  !NRT
     Orbit = trim(AMSRFILE(L1R_dir_len+24:L1R_dir_len+24))
  elseif(AMSReOPL%L1Rtype.eq.2) then  !Historical
     Orbit = trim(AMSRFILE(L1R_dir_len+20:L1R_dir_len+20))
  endif

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
  !deallocate(SCNANG_L1B)
  !deallocate(ANTSCN_L1B)
  deallocate(FLAG_L1R)
  !deallocate(TBHFLAG_L1B)
  deallocate(ARFS_LAT)
  deallocate(ARFS_LON)

end subroutine AMSR_L1R_SUBSET
