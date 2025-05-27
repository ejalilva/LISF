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
! SUBROUTINE: ARFSSMRETRIEVAL
!
! REVISION HISTORY:
!  22 Feb 2022: P.W.LIU; Initial implemetation
!  22 Feb 2022: Yonghwan Kwon; modified for LDT
!  10 Feb 2023: Eric Kemp, modified to output retrievals in netCDF.
!  21 Feb 2023: Eric Kemp, added third LIS time level.
!
! DESCRIPTION: RETRIEVE AMSR SM FOR ARFS
! INPUT : AMSR - L1R Brightness Temperature
! OUTPUT: AMSRTB_ARFSGRIDE_ddmmyyy.dat
! NOTES : Inverse Distance Squared with 0.4 deg serching window
!-------------------------------------------------------------------------

subroutine ARFSSMRETRIEVAL_AMSR(AMSRFILE, &
     TS_bfresample_01, TS_bfresample_02, TS_bfresample_03, &
     ARFS_SNOW, DOY, UTChr, firsttime, secondtime, thirdtime)

   !USE HDF5
    use esmf
    USE VARIABLES
    USE DATADOMAIN
    USE FUNCTIONS
    USE netcdf
    USE invdist_temp2amsr ! E.J: this is for resampling to smap 33 km grid, do we need to do it for AMSR?
    USE varsio_m_amsr
    USE algo_hpol_m
    use LDT_AMSR_ARFSSM_netcdfMod, only: LDT_AMSR_ARFSSM_write_netcdf ! E.J: check what is this for
    use LDT_logMod, only: LDT_logunit
    USE LDT_amsr_oplMod  

    IMPLICIT NONE
! !ARGUMENTS:
    CHARACTER (len=100)          :: AMSRFILE                             
    REAL*4, DIMENSION(2560,1920), intent(in) :: TS_bfresample_01, &
         TS_bfresample_02, TS_bfresample_03
    REAL*4, DIMENSION(2560,1920) :: ARFS_SNOW, UTChr                     
    INTEGER                      :: DOY                                  
    type(ESMF_Time), intent(in) :: firsttime
    type(ESMF_Time), intent(in) :: secondtime
    type(ESMF_Time), intent(in) :: thirdtime
!EOP 
    INTEGER :: i, j, nrow, mcol          
    CHARACTER (len=100) :: fname_TAU    
    CHARACTER (len=5) :: DOY_chr
    REAL*4 :: C, K, sm_retrieval, tau_return
    REAL*4, DIMENSION(2560,1920) :: ARFS_TB, ARFS_TB_37V
    REAL*4, DIMENSION(2560,1920) :: ARFS_TAU, ARFS_CLAY, ARFS_BD, ARFS_OMEGA, ARFS_H
    INTEGER*1, DIMENSION(2560,1920) :: ARFS_LC, ARFS_SM_FLAG
    INTEGER*1 :: retrieval_flag
    REAL*4, DIMENSION(2560,1920) :: ARFS_TS_01, ARFS_TS_02, ARFS_TS_03, ARFS_SM
    REAL*8 ,DIMENSION(:), ALLOCATABLE :: ARFS_FINE_LAT, ARFS_FINE_LON
    REAL*8 ,DIMENSION(:), ALLOCATABLE :: ARFS_LAT, ARFS_LON
    REAL :: T1, T2
    
    INTEGER*4 :: ios, NX, NY
    INTEGER*4 :: ncid, nid, tsoil01id

    character (len=100) :: retrieval_fname
    integer             :: L1R_dir_len,L1R_fname_len
    real                :: utc_check

    ! EMK
    character(8) :: yyyymmdd
    character(6) :: hhmmss
    real :: deltasec, wgt
    integer :: firstUTCyr, firstUTCmo, firstUTCdy, firstUTChr
    integer :: secondUTCyr, secondUTCmo, secondUTCdy, secondUTChr
    integer :: thirdUTCyr, thirdUTCmo, thirdUTCdy, thirdUTChr

    real :: TS_A, TS_B

    integer :: count_valid_points = 0
    integer :: count_tbh_invalid = 0
    integer :: count_ts_invalid = 0
    integer :: count_snow_invalid = 0
    integer :: count_bd_invalid = 0
    integer :: count_lc_invalid = 0
    integer :: count_utc_invalid = 0

    nrow=2560
    mcol=1920

    ! RESAMPLE TEFF TO 33 KM on ARFS GRID
    write (LDT_logunit,*) '[INFO] Resampling effective soil temperature'
    CALL ARFS_GEO
    ALLOCATE(ARFS_LAT(arfs_nrow_lat),ARFS_LON(arfs_mcol_lon))
    ARFS_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_space)
    ARFS_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_space)
    CALL ARFS_3KM_GEO
    ALLOCATE(ARFS_FINE_LAT(arfs_nrow_lat_3km),ARFS_FINE_LON(arfs_mcol_lon_3km))
    ARFS_FINE_LAT = LAT(arfs_geo_lat_lo,arfs_geo_lat_up,-arfs_lat_3km_space)
    ARFS_FINE_LON = LON(arfs_geo_lon_lf,arfs_geo_lon_rt,arfs_lon_3km_space)
    CALL RESAMPLETEMP(TS_bfresample_01,ARFS_LAT,ARFS_LON,ARFS_FINE_LAT,ARFS_FINE_LON,ARFS_TS_01)
    CALL RESAMPLETEMP(TS_bfresample_02,ARFS_LAT,ARFS_LON,ARFS_FINE_LAT,ARFS_FINE_LON,ARFS_TS_02)
    CALL RESAMPLETEMP(TS_bfresample_03,ARFS_LAT,ARFS_LON,ARFS_FINE_LAT,ARFS_FINE_LON,ARFS_TS_03)
    ! IF EVENTUALLY THE RESAMPLING DOES NOT CHANGE TEFF MUCH WE COULD SIMPLELY USE ARFS_TS=TS_bfresample
    ! UP TO HERE TAKES 38 SECS
    write (LDT_logunit,*) '[INFO] Finished resampling effective soil temperature'

    ! get RESAMPLED TB
    ARFS_TB = AMSReOPL%ARFS_TB_10H
    ARFS_TB_37V = AMSReOPL%ARFS_TB_36V ! added for calculation of TS based on 37GHz band

    ! LOAD TAU ------------------------------------------------------------
    write(DOY_chr,"(I0.3)") DOY
    fname_TAU = trim(AMSReOPL%TAUdir)//"/tau_cmg_arfs_"//trim(DOY_chr)//".dat"
    write (LDT_logunit,*) '[INFO] Reading TAU from ', trim(fname_TAU)
    OPEN(UNIT=1,FILE=fname_TAU,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD',convert='little_endian')
    READ(1, rec=1) ARFS_TAU
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading TAU'

    ! LOAD OMEGA
    write (LDT_logunit,*) '[INFO] Reading OMEGA from ', trim(AMSReOPL%OMEGAfile)
    OPEN(UNIT=1,FILE=AMSReOPL%OMEGAfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD', convert='little_endian')
    READ(1, rec=1) ARFS_OMEGA
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading OMEGA'

    ! LOAD SOIL (BD and clayFraction)
    write (LDT_logunit,*) '[INFO] Reading soil bulk density from ', trim(AMSReOPL%BDfile)
    OPEN(UNIT=1,FILE=AMSReOPL%BDfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD', convert='little_endian') !Bulk Density
    READ(1, rec=1) ARFS_BD
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading soil bulk density'

    ! E.J: Debugging ********
    write(LDT_logunit,*) '[INFO] PARAMETER STATISTICS:'
    write(LDT_logunit,*) '  Bulk Density - min:', MINVAL(ARFS_BD), 'max:', MAXVAL(ARFS_BD)
    write(LDT_logunit,*) '  Calculated upperbound range:', MINVAL(1-ARFS_BD/2.65), 'to', MAXVAL(1-ARFS_BD/2.65)

    
    write (LDT_logunit,*) '[INFO] Reading soil clay fraction from ', trim(AMSReOPL%CLAYfile)
    OPEN(UNIT=1,FILE=AMSReOPL%CLAYfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD',convert='little_endian') !Clay Fraction
    READ(1, rec=1) ARFS_CLAY
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading soil clay fraction'

    ! LOAD ROUGHNESS
    write (LDT_logunit,*) '[INFO] Reading roughness from ', trim(AMSReOPL%Hfile)
    OPEN(UNIT=1,FILE=AMSReOPL%Hfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*nrow*mcol,STATUS='OLD',convert='little_endian') !roughness
    READ(1, rec=1) ARFS_H
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading roughness'

    ! LOAD LANDCOVER
    write (LDT_logunit,*) '[INFO] Reading landcover from ', trim(AMSReOPL%LCfile)
    OPEN(UNIT=1,FILE=AMSReOPL%LCfile,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1*nrow*mcol,STATUS='OLD',convert='little_endian')
    READ(1, rec=1) ARFS_LC
    CLOSE(1)
    write (LDT_logunit,*) '[INFO] Finished reading landcover'

    !generate soil moisture retrievals
    write (LDT_logunit,*) '[INFO] Generating soil moisture retrievals'
    ARFS_SM=-9999
    ARFS_SM_FLAG=-1

    call ESMF_TimeGet(firsttime, yy=firstUTCyr, mm=firstUTCmo, dd=firstUTCdy, &
         h=firstUTChr)
    call ESMF_TimeGet(secondtime, yy=secondUTCyr, mm=secondUTCmo, &
         dd=secondUTCdy, &
         h=secondUTChr)
    call ESMF_TimeGet(thirdtime, yy=thirdUTCyr, mm=thirdUTCmo, dd=thirdUTCdy, &
         h=thirdUTChr)
    
    ! E.J: logging for debugging purpose
    write(LDT_logunit,*) '[DEBUG] Time interpolation:'
    write(LDT_logunit,*) '  firstUTCyr, firstUTCmo, firstUTCdy, firstUTChr:', &
         firstUTCyr, firstUTCmo, firstUTCdy, firstUTChr
    write(LDT_logunit,*) '  secondUTCyr, secondUTCmo, secondUTCdy, secondUTChr:', &
         secondUTCyr, secondUTCmo, secondUTCdy, secondUTChr
    write(LDT_logunit,*) '  thirdUTCyr, thirdUTCmo, thirdUTCdy, thirdUTChr:', &
         thirdUTCyr, thirdUTCmo, thirdUTCdy, thirdUTChr

    write(LDT_logunit,*) '  ARFS_TS_01 stats (min, max, count<0):', &
         MINVAL(ARFS_TS_01, MASK=(ARFS_TS_01 > -9990)), &
         MAXVAL(ARFS_TS_01, MASK=(ARFS_TS_01 > -9990)), &
         COUNT(ARFS_TS_01 <= 0)
    write(LDT_logunit,*) '  ARFS_TS_02 stats (min, max, count<0):', &
         MINVAL(ARFS_TS_02, MASK=(ARFS_TS_02 > -9990)), &
         MAXVAL(ARFS_TS_02, MASK=(ARFS_TS_02 > -9990)), &
         COUNT(ARFS_TS_02 <= 0)
    write(LDT_logunit,*) '  ARFS_TS_03 stats (min, max, count<0):', &
         MINVAL(ARFS_TS_03, MASK=(ARFS_TS_03 > -9990)), &
         MAXVAL(ARFS_TS_03, MASK=(ARFS_TS_03 > -9990)), &
         COUNT(ARFS_TS_03 <= 0)
     
    DO j=1,mcol !COL LAT
       DO i=1,nrow !ROW LON

          tbh = ARFS_TB(i,j) 

          if (UTChr(i,j) < 0) cycle

          if (UTChr(i,j) == firstUTChr) then
             TS_A = ARFS_TS_01(i,j)
             TS_B = ARFS_TS_02(i,j)
             wgt = 1
          else if (UTChr(i,j) > firstUTChr .and. &
               UTChr(i,j) < secondUTChr) then
             TS_A = ARFS_TS_01(i,j)
             TS_B = ARFS_TS_02(i,j)
             deltasec = ( UTChr(i,j) - firstUTChr ) * 3600
             wgt = (10800. - deltasec) / 10800.
          else if (UTChr(i,j) > firstUTChr .and. &
               firstUTChr == 21 .and. secondUTChr == 0) then
             TS_A = ARFS_TS_01(i,j)
             TS_B = ARFS_TS_02(i,j)
             deltasec = ( UTChr(i,j) - firstUTChr ) * 3600
             wgt = (10800. - deltasec) / 10800.
          else if (UTChr(i,j) == secondUTChr) then
             TS_A = ARFS_TS_02(i,j)
             TS_B = ARFS_TS_03(i,j)
             wgt = 1
          else
             TS_A = ARFS_TS_02(i,j)
             TS_B = ARFS_TS_03(i,j)
             deltasec = ( UTChr(i,j) - secondUTChr ) * 3600
             wgt = (10800. - deltasec) / 10800.
          end if
          if (TS_A > 0 .and. TS_B > 0) then
             TS = ((wgt)*TS_A) + ((1. - wgt)*TS_B) ! E.J: This has just being used to test if TS is above 0 degree
          else
             cycle
          end if
          
          ! ========
          ! E.J: for debugging purpose just use the TS_01 to avoid using scan time data that is buggy and see if SM retrieval is working
          !TS = ARFS_TS_01(i,j)
          ! E.J: using the 37 GHz band TS 
          TS =1.11*ARFS_TB_37V(i,j)-15.2
          !==========
          
          !IF (tbh.GT.0.0.AND.Ts.GT.0.AND.ARFS_SNOW(i,j).LE.AMSReOPL%SD_thold.AND.ARFS_BD(i,j).NE.-9999.AND.ARFS_LC(i,j).NE.0.AND.&
           ! UTChr(i,j).GE.0) THEN
          IF (tbh.LE.0.0) THEN
             count_tbh_invalid = count_tbh_invalid + 1
          ELSEIF (Ts.LE.0) THEN
             count_ts_invalid = count_ts_invalid + 1
          ELSEIF (ARFS_SNOW(i,j).GT.AMSReOPL%SD_thold) THEN
             count_snow_invalid = count_snow_invalid + 1
          ELSEIF (ARFS_BD(i,j).EQ.-9999) THEN
             count_bd_invalid = count_bd_invalid + 1
          ELSEIF (ARFS_LC(i,j).EQ.0) THEN
             count_lc_invalid = count_lc_invalid + 1
          ELSEIF (UTChr(i,j).LT.0) THEN
             count_utc_invalid = count_utc_invalid + 1
          ELSE
             count_valid_points = count_valid_points + 1
             bulkdensity = ARFS_BD(i,j)
             clay = ARFS_CLAY(i,j)
             tau = ARFS_TAU(i,j)*1.3353 ! changing to account for conversion from L-band to xband to account for frequency difference multiplied by Sec(55)/sec(40) ~ 1.3353
             omega = ARFS_OMEGA(i,j) ! changing to accout for conversion from L-band to xband 
             h = ARFS_H(i,j)
             topigbptype = ARFS_LC(i,j)

             CALL algo_hpol(real(i),real(j),sm_retrieval, tau_return, retrieval_flag)
             ARFS_SM(i,j)=sm_retrieval
             ARFS_SM_FLAG(i,j)=retrieval_flag

          END IF
       END DO !ii=1,nrow !ROW LON
    END DO !jj=1,mcol !COL LAT

    !write soil moisture retrieval outputs
    L1R_dir_len = len_trim(AMSReOPL%L1Rdir)
    L1R_fname_len = len_trim(AMSRFILE)

    ! TODO: make H and V automatic based on which pol is used in retrieval (also LIS teff or 37teff should be automatically reflected in the naming.
    
    if(AMSReOPL%L1Rtype.eq.1) then  !NRT
       retrieval_fname = trim(AMSReOPL%SMoutdir)//"/"//"ARFS_SM_H_"//& 
                         trim(AMSRFILE(L1R_dir_len+9:L1R_fname_len-3))//".nc"
       yyyymmdd = trim(AMSRFILE(L1R_dir_len+9:L1R_dir_len+15))
       hhmmss = trim(AMSRFILE(L1R_dir_len+16:L1R_dir_len+19))
    elseif(AMSReOPL%L1Rtype.eq.2) then  !Historical
       retrieval_fname = trim(AMSReOPL%SMoutdir)//"/"//"ARFS_SM_H_"//&
            trim(AMSRFILE(L1R_dir_len+9:L1R_fname_len-3))//".nc"
       yyyymmdd = trim(AMSRFILE(L1R_dir_len+9:L1R_dir_len+15))
       hhmmss = trim(AMSRFILE(L1R_dir_len+16:L1R_dir_len+19))
       write (LDT_logunit,*) 'yyyymmhh: ', trim(yyyymmdd), ', hhmmss: ', trim(hhmmss)

    endif

    write (LDT_logunit,*) '[INFO] Writing soil moisture retrieval file ', trim(retrieval_fname)

    ! NOTE: nrow is actually number of columns, mcol is actually number of
    ! rows
    call LDT_AMSR_ARFSSM_write_netcdf(nrow, mcol, arfs_sm, retrieval_fname, &
         yyyymmdd, hhmmss)
    write (LDT_logunit,*) '[INFO] Successfully wrote soil moisture retrieval file ', trim(retrieval_fname)
    write (LDT_logunit,*) '[INFO] Finished generating soil moisture retrievals'
    write(LDT_logunit,*) '[DEBUG] Filtering statistics:'
    write(LDT_logunit,*) '  Valid points processed:', count_valid_points
    write(LDT_logunit,*) '  Points with invalid tbh:', count_tbh_invalid
    write(LDT_logunit,*) '  Points with invalid Ts:', count_ts_invalid
    write(LDT_logunit,*) '  Points with snow above threshold:', count_snow_invalid
    write(LDT_logunit,*) '  Points with invalid bulk density:', count_bd_invalid
    write(LDT_logunit,*) '  Points with invalid land cover:', count_lc_invalid
    write(LDT_logunit,*) '  Points with invalid UTC hour:', count_utc_invalid

 end subroutine ARFSSMRETRIEVAL_AMSR
