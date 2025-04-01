!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE algo_hpol_m
        IMPLICIT NONE
 
      CONTAINS
       FUNCTION algo_hpol_function(X) RESULT(simtbh)
        USE varsio_m_amsr
        USE mironov_m

        IMPLICIT NONE

        REAL*4    :: roh, rov, rsh, rsv, exptauh, exptauv, Ah, Av, X
        COMPLEX(4) :: er_r, c_er, er
        REAL*4    :: algo_hpol_output, simtbv, simtbh

        CALL mironov (freq,X,clay,er_r) !freq: frequency in GHZ

        c_er = er_r
        er = SQRT (c_er - SIN (inc*d2r)**2)
        roh = ABS ( (COS (inc*d2r) - er) / (COS (inc*d2r) + er) )**2
        rov = ABS ( (c_er * COS (inc*d2r) - er) / (c_er * COS (inc*d2r) + er) )**2
        rsh = ( (1-Q) * roh + Q * rov) * EXP (-h * COS (inc*d2r) * COS (inc*d2r))
        rsv = ( (1-Q) * rov + Q * roh) * EXP (-h * COS (inc*d2r) * COS (inc*d2r))
        exptauh = EXP (-tau)
        exptauv = EXP (-tau)
        Ah = Ts * (1 - omega) * (1 - exptauh)
        Av = Ts * (1 - omega) * (1 - exptauv)
        simtbh = Ts * (1 - rsh) * exptauh + Ah * (1 + rsh * exptauh)
        simtbv = Ts * (1 - rsv) * exptauv + Av * (1 + rsv * exptauv)
        algo_hpol_output = simtbh

      END FUNCTION algo_hpol_function

        SUBROUTINE algo_hpol (ii,jj,x1,x2,exitstate)
          USE varsio_m_amsr
          IMPLICIT NONE
 
          REAL(4), INTENT(IN)                      :: ii, jj
          REAL(4), INTENT(OUT)                     :: x1, x2
          INTEGER(1), INTENT(OUT)                  :: exitstate
          REAL(4)                                  :: NEDT = 2.0
          REAL(4)                                  :: lowerbound = 0.02
          REAL(4)                                  :: upperbound
          REAL(4)                                  :: x
          REAL(4)                                  :: incvsm
          INTEGER(4)                               :: numvsm
          INTEGER(4)                               :: hh, opt
          REAL(4), DIMENSION(:), ALLOCATABLE, SAVE :: vsmvec, tbhvec
          incvsm = 0.01
          upperbound = 1 - bulkdensity/2.65
          numvsm = FLOOR ((upperbound - lowerbound)/incvsm)
          ALLOCATE (vsmvec(numvsm))
          ALLOCATE (tbhvec(numvsm))
          DO hh = 1,numvsm
             vsmvec(hh) = lowerbound + (numvsm-1)*incvsm - (hh-1)*incvsm
             tbhvec(hh) = algo_hpol_function(vsmvec(hh))
          ENDDO
          IF (tbh >= tbhvec(1) - NEDT .AND. tbh <= tbhvec(numvsm) + NEDT) THEN
              IF (tbh < tbhvec(1)) THEN ! assigning tbh of residual soil moisture if tbh is smaller than smallest tbh
                  tbh = tbhvec(1)
              ENDIF
              IF (tbh > tbhvec(numvsm)) THEN ! assigning tbh of saturated soil moisture if tbh is higher than highest tbh
                  tbh = tbhvec(numvsm)
              ENDIF
               IF (tbhvec(numvsm) - tbhvec(1) > NEDT) THEN
                   opt=MINLOC(ABS(tbh-tbhvec),1)
                   x = vsmvec(opt)
                   exitstate = 0
               ELSE
                   x = FillValue_float32
                   exitstate = 1
               ENDIF
          ELSEIF (tbh > tbhvec(numvsm) + NEDT) THEN
              IF (topigbptype >= 1 .AND. topigbptype <= 5) THEN
                  x = upperbound
              ELSE
                  x = lowerbound
              ENDIF
              exitstate = 1
          ELSEIF (tbh < tbhvec(1) - NEDT) THEN
              x = upperbound
              exitstate = 1
          ELSE
              x = FillValue_float32
              exitstate = 1
          ENDIF
          x1 = x
          x2 = tau
          DEALLOCATE (vsmvec)
          DEALLOCATE (tbhvec)
 
        END SUBROUTINE algo_hpol
 
      END MODULE algo_hpol_m
