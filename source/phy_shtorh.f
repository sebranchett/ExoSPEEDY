 
      SUBROUTINE SHTORH (IMODE,NGP,TA,PS,SIG,QA,RH,QSAT)
C--
C--   SUBROUTINE SHTORH (IMODE,NGP,TA,PS,SIG,QA,RH,QSAT)
C--
C--   Purpose: compute saturation specific humidity and 
C--            relative hum. from specific hum. (or viceversa)
C--   Input:   IMODE  : mode of operation
C--            NGP    : no. of grid-points
C--            TA     : abs. temperature
C--            PS     : normalized pressure   (=  p/1000_hPa) [if SIG < 0]
C--                   : normalized sfc. pres. (= ps/1000_hPa) [if SIG > 0]
C--            SIG    : sigma level
C--            QA     : specific humidity in g/kg [if IMODE > 0]
C--            RH     : relative humidity         [if IMODE < 0]
C--            QSAT   : saturation spec. hum. in g/kg 
C--   Output:  RH     : relative humidity         [if IMODE > 0] 
C--            QA     : specific humidity in g/kg [if IMODE < 0]
C--        
C--IO h planetparam.h, com_planet.h
C--IO sx E0=  6.108E-3 - August-Roche-Magnus factor ARMFAC
C--IO sx C1= 17.269 - August-Roche-Magnus constant above freezing ARMC1
C--IO sx C2= 21.875 - August-Roche-Magnus constant below freezing ARMC2
C--IO sx 622. - ratio water vapor to dry air with unit correction
C--IO sx 0.378 = 1-.622
C--IO sx T0=273.16 - temperature in degrees Kelvin 0C
C--IO sx T1= 35.86 - August-Roche-Magnus temperature above freezing ARMT1
C--IO sx T2=  7.66 - August-Roche-Magnus temperature below freezing ARMT2
      include "planetparam.h"
      include "com_planet.h"

      REAL TA(NGP), PS(*), QA(NGP), RH(NGP), QSAT(NGP)
C
C---  1. Compute Qsat (g/kg) from T (degK) and normalized pres. P (= p/1000_hPa)
C        If SIG > 0, P = Ps * sigma, otherwise P = Ps(1) = const. 
C
      E0= ARMFAC
      C1= ARMC1
      C2= ARMC2
      T0= FRWTR2
      T1= ARMT1
      T2= ARMT2
      
      DO 110 J=1,NGP
        IF (TA(J).GE.T0) THEN
          QSAT(J)=E0*EXP(C1*(TA(J)-T0)/(TA(J)-T1))
        ELSE
          QSAT(J)=E0*EXP(C2*(TA(J)-T0)/(TA(J)-T2))
        ENDIF
  110 CONTINUE
C
C     Eq. (41) of:
C     Bolton, D. (1980) The Computation of Equivalent Potential
C     Temperature. Monthly Weather Review, 108, 1046-1053.
C     https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2
C     see also: https://cran.r-project.org/web/packages/humidity/vignettes/
C               humidity-measures.html
      ONEMIN = 1. - WTRAIR
      IF (SIG.LE.0.0) THEN
        DO 120 J=1,NGP
          QSAT(J)=1000*WTRAIR*QSAT(J)/(PS(1)-ONEMIN*QSAT(J))
  120   CONTINUE
      ELSE
        DO 130 J=1,NGP
          QSAT(J)=1000*WTRAIR*QSAT(J)/(SIG*PS(J)-ONEMIN*QSAT(J))
  130   CONTINUE
      ENDIF
C
C---  2. Compute rel.hum. RH=Q/Qsat (IMODE>0), or Q=RH*Qsat (IMODE<0)
C
      IF (IMODE.GT.0) THEN
        DO 210 J=1,NGP
          RH(J)=QA(J)/QSAT(J)
  210   CONTINUE
      ELSE IF (IMODE.LT.0) THEN
        DO 220 J=1,NGP
          QA(J)=RH(J)*QSAT(J)
  220   CONTINUE
      ENDIF
C				
      RETURN
      END

      SUBROUTINE ZMEDDY (NLON,NLAT,FF,ZM,EDDY)
C
C *** Decompose a field into zonal-mean and eddy component
C
      REAL FF(NLON,NLAT), ZM(NLAT), EDDY(NLON,NLAT)
C
      RNLON=1./NLON
C
      DO 130 J=1,NLAT
C
        ZM(J)=0.
        DO 110 I=1,NLON
          ZM(J)=ZM(J)+FF(I,J)
 110    CONTINUE
        ZM(J)=ZM(J)*RNLON
C
        DO 120 I=1,NLON
          EDDY(I,J)=FF(I,J)-ZM(J)
 120    CONTINUE
C
 130  CONTINUE
C
C--
      RETURN
      END




