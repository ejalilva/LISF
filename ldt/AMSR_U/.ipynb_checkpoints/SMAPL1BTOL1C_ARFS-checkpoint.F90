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
! SUBROUTINE: SMAPL1BRESAMPLE
!
! REVISION HISTORY:
!  22 Oct 2021: P.W.LIU; Initial implemetation
!  18 Dec 2021: Yonghwan Kwon; modified for LDT
!
! DESCRIPTION: RESAMPLE SMAPL1B TB TO AIR FORCE GRID
! INPUT : SMAP - L1B Brightness Temperature               
! OUTPUT: SMAPTB_ARFSGRIDE_ddmmyyy.dat
! NOTES : Inverse Distance Squared with 0.4 deg serching window
!-------------------------------------------------------------------------

subroutine AMSR_L1R_RESAMPLE(AMSRFILE,L1R_dir,Orbit,ARFS_TIME,rc)

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
    character (len=20)  :: variable_name(13)
    character (len=100) :: resample_filename(13)
    character (len=1)   :: Orbit
    integer             :: var_i
    integer             :: L1R_dir_len,L1R_fname_len
    integer             :: ierr
    integer             :: rc

    REAL*8,DIMENSION(:,:),ALLOCATABLE :: TIME_L1R
    REAL*4,DIMENSION(:,:),ALLOCATABLE :: TB_10H, TB_10V, TB_18H, TB_18V, TB_23H, TB_23V, TB_36H, TB_36V, TB_89H, TB_89V, LAND_WATER_FRAC
    !REAL*4,DIMENSION(:,:),ALLOCATABLE :: NETD_V_L1B, NETD_H_L1B, LAT_L1B, LON_L1B, SCNANG_L1B
    !REAL*4,DIMENSION(:),ALLOCATABLE :: ANTSCN_L1B
    INTEGER*4,DIMENSION(:,:),ALLOCATABLE :: FLAG_L1R


    REAL*8,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
    REAL*8,DIMENSION(2560,1920) :: ARFS_TIME
    !INTEGER*4,DIMENSION(2560,1920) :: ARFS_SAMPLE_V, ARFS_SAMPLE_H 
    REAL*4,DIMENSION(2560,1920) :: ARFS_TB_10V, ARFS_TB_10H, ARFS_TB_18H, ARFS_TB_18V,ARFS_TB_23H, ARFS_TB_23V, ARFS_TB_36H, ARFS_TB_36V, ARFS_TB_89H, ARFS_TB_89V, ARFS_LAND_WATER_FRAC
    INTEGER*4,DIMENSION(2560,1920),ALLOCATABLE :: ARFS_FLAG_L1R

    !REAL*4,DIMENSION(2560,1920) :: ARFS_SURWAT_V, ARFS_SURWAT_H, ARFS_WTV, ARFS_WTH

    REAL :: T1, T2
    CALL ARFS_GEO
    ALLOCATE(ARFS_LAT(arfs_nrow_lat),ARFS_LON(arfs_mcol_lon))
    ARFS_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_space)
    ARFS_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_space)
    
    !Input (Path/filename, datatypes, lat, lon) Return(DATA,LAT,LON,length of row and col) 
    ! READ SMAP_L1B DATA FROM HDF5

    CALL GetAMSR_L1R(AMSRFILE,TIME_L1R &
            TB_10H, TB_10V, TB_18H, TB_18V, &
            TB_23H, TB_23V, TB_36H, TB_36V, &
            TB_89H, TB_89V, &
            LAT_L1R, LON_L1R, &
            FLAG_L1R, LAND_WATER_FRAC, &
            nrow, mcol, ierr)

    !Input (DATA,LAT,LON,length of row and col); Return(TB in ARFS GRID)
    CALL L1BTB2ARFS_INVDIS(TIME_L1B, TBV_COR_L1B, TBH_COR_L1B, TBV_L1B, TBH_L1B, SURWAT_V_L1B, SURWAT_H_L1B, &
                           NETD_V_L1B, NETD_H_L1B, LAT_L1B, LON_L1B, TBVFLAG_L1B, TBHFLAG_L1B, ANTSCN_L1B, SCNANG_L1B, nrow, mcol, &
                           ARFS_LAT, ARFS_LON, ARFS_TIME, ARFS_COR_TBV, ARFS_COR_TBH, ARFS_TBV, ARFS_TBH, ARFS_NEDTV, ARFS_NEDTH, &
                           ARFS_SURWAT_V, ARFS_SURWAT_H, ARFS_WTV, ARFS_WTH, ARFS_SAMPLE_V, ARFS_SAMPLE_H)

    SMAPeOPL%ARFS_TBV_COR = ARFS_COR_TBV

    variable_name(1)  = 'ARFS_TIME'
    variable_name(2)  = 'ARFS_TBV_COR'
    variable_name(3)  = 'ARFS_TBH_COR'
    variable_name(4)  = 'ARFS_TBV'
    variable_name(5)  = 'ARFS_TBH'
    variable_name(6)  = 'ARFS_NETDV'
    variable_name(7)  = 'ARFS_NETDH'
    variable_name(8)  = 'ARFS_SURWATV'
    variable_name(9)  = 'ARFS_SURWATH'
    variable_name(10) = 'ARFS_WTV'
    variable_name(11) = 'ARFS_WTH'
    variable_name(12) = 'ARFS_SAMPV'
    variable_name(13) = 'ARFS_SAMPH'

    L1B_dir_len = len_trim(L1B_dir)
    L1B_fname_len = len_trim(SMAPFILE)

    if(SMAPeOPL%L1Btype.eq.1) then  !NRT
       Orbit = trim(SMAPFILE(L1B_dir_len+24:L1B_dir_len+24))
    elseif(SMAPeOPL%L1Btype.eq.2) then  !Historical
       Orbit = trim(SMAPFILE(L1B_dir_len+20:L1B_dir_len+20))
    endif

    if(SMAPeOPL%L1BresampWriteOpt.eq.1) then
       if(SMAPeOPL%L1Btype.eq.1) then  !NRT
          do var_i=1,13
             resample_filename(var_i) = trim(SMAPeOPL%L1Bresampledir_02)//"/"//trim(variable_name(var_i))//"_"//&
                                        trim(SMAPFILE(L1B_dir_len+18:L1B_fname_len-3))//".dat"
          enddo
       elseif(SMAPeOPL%L1Btype.eq.2) then  !Historical
          do var_i=1,13
             resample_filename(var_i) = trim(SMAPeOPL%L1Bresampledir_02)//"/"//trim(variable_name(var_i))//"_"//&
                                        trim(SMAPFILE(L1B_dir_len+14:L1B_fname_len-3))//".dat"
          enddo
       endif

       OPEN(UNIT=151, FILE=resample_filename(1),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TIME
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(2),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_COR_TBV
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(3),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_COR_TBH
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(4),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TBV
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(5),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_TBH
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(6),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_NEDTV
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(7),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_NEDTH
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(8),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_SURWAT_V
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(9),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_SURWAT_H
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(10),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_WTV
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(11),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_WTH
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(12),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_SAMPLE_V
       CLOSE(151)
       OPEN(UNIT=151, FILE=resample_filename(13),FORM='UNFORMATTED',ACCESS='DIRECT', RECL=arfs_nrow_lat*arfs_mcol_lon*4)
       WRITE(UNIT=151, REC = 1) ARFS_SAMPLE_H
       CLOSE(151)
    endif

end subroutine SMAPL1BRESAMPLE
