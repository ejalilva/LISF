!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE varsio_m_amsr
       IMPLICIT NONE
       REAL*8, PARAMETER :: RE_KM = 6371.228, PI = acos(-1.0), d2r = PI/180.0
       REAL*4, PARAMETER :: FillValue_float32=-9999, Q=0, freq=10.65, inc=55!53.8
       REAL*4 :: bulkdensity, tbh, tau, clay, h, Ts, omega
       INTEGER*1 :: topigbptype 
END MODULE varsio_m_amsr
