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
  USE TOOLSUBS
  USE invdist_l1r2arfs
  USE LDT_logMod
  USE LDT_amsr_oplMod

  IMPLICIT NONE

  INTEGER :: i, j, nrow, mcol
  CHARACTER (len=100) :: AMSRFILE
  character (len=100) :: L1R_dir
  character (len=20)  :: variable_name(12)
  character (len=100) :: resample_filename(12)
  character (len=1)   :: Orbit ! E.J: orbit is one of the outputs extracted from the filename (D: Descending & A: Ascending)
  integer             :: var_i
  integer             :: L1R_dir_len,L1R_fname_len
  integer :: ierr
  integer :: rc

  REAL*8,DIMENSION(:,:),ALLOCATABLE :: TIME_L1R
  REAL*4,DIMENSION(:,:),ALLOCATABLE :: TB_10H, TB_10V, TB_18H, TB_18V, TB_23H, TB_23V, TB_36H, TB_36V, TB_89H, TB_89V, LAND_WATER_FRAC
  REAL*4,DIMENSION(:,:),ALLOCATABLE :: LAT_L1R, LON_L1R, LAT89, LON89
  INTEGER*4,DIMENSION(:,:),ALLOCATABLE :: RFI_FLAG
  INTEGER*4,DIMENSION(:,:),ALLOCATABLE :: SNOW, PRECIP
  INTEGER :: nrow89,ncol89 ! m89 & n89 are for the 89GHz band for which lat and lon are provided

  REAL*8,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
  INTEGER*4,DIMENSION(2560,1920) :: ARFS_SAMPLE_V, ARFS_SAMPLE_H 
  REAL*8,DIMENSION(2560,1920) :: ARFS_TIME
  REAL*4,DIMENSION(2560,1920) :: ARFS_TB_10V, ARFS_TB_10H, ARFS_TB_18H, ARFS_TB_18V,ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, ARFS_TB_89H, ARFS_TB_89V, ARFS_LAND_WATER_FRAC ! This is instead of ARFS_COR_TBV
  !INTEGER*4,DIMENSION(2560,1920),ALLOCATABLE :: ARFS_RFI_FLAG

  REAL :: T1, T2

  rc = 0

  CALL ARFS_GEO
  ALLOCATE(ARFS_LAT(arfs_nrow_lat),ARFS_LON(arfs_mcol_lon)) 
  ARFS_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_space)
  ARFS_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_space)

  CALL get_amsr_l1r (AMSRFILE,TIME_L1R, &
          TB_10V, TB_10H, TB_18V, TB_18H, &
          TB_23V, TB_23H, TB_36V, TB_36H, &
          TB_89V, TB_89H, &
          LAT_L1R, LON_L1R, LAT89, LON89, &
          LAND_WATER_FRAC, SNOW, PRECIP, &
          RFI_FLAG, &
          nrow, mcol, nrow89,ncol89, ierr)
  ! TODO check all of these and find where they are called, nrow89,ncol89 perhaps in the invdist script the input to invdist should be modified
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
          SNOW, PRECIP, &
          LAT_L1R, LON_L1R, nrow, mcol, &
          ARFS_LAT, ARFS_LON, ARFS_TIME, ARFS_LAND_WATER_FRAC, &
          ARFS_TB_10H, ARFS_TB_10V, ARFS_TB_18H, ARFS_TB_18V, &
          ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, &
          ARFS_TB_89H, ARFS_TB_89V, ARFS_SAMPLE_V, ARFS_SAMPLE_H)
  
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
    if(AMSReOPL%L1RresampWriteOpt.eq.1) then
       if(AMSReOPL%L1Rtype.eq.1) then  !NRT
          do var_i=1,12
             resample_filename(var_i) = trim(AMSReOPL%L1Rresampledir_02)//"/"//trim(variable_name(var_i))//"_"//& ! EJ: Where L1Bresampledir_02 is being set
                                        trim(AMSRFILE(L1R_dir_len+18:L1R_fname_len-3))//".dat"
          enddo
       elseif(AMSReOPL%L1Rtype.eq.2) then  !Historical
          do var_i=1,12
             resample_filename(var_i) = trim(AMSReOPL%L1Rresampledir_02)//"/"//trim(variable_name(var_i))//"_"//&
                                        trim(AMSRFILE(L1R_dir_len+14:L1R_fname_len-3))//".dat"
          enddo
       endif

       OPEN(UNIT=151, FILE=resample_filename(1),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TIME
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(2),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_10H
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(3),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_10V
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(4),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_18H
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(5),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_18V
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(6),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_23H
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(7),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_23V
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(8),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_36H
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(9),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_36V
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(10),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_89H
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(11),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TB_89V
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(12),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_LAND_WATER_FRAC
       CLOSE(151)
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

end subroutine AMSR_L1R_RESAMPLE
