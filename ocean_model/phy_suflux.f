
      SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
     &                   PHI0,TSEA,SSRD,SLRD,
     &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
     &                   TSFC,TSKIN,U0,V0,T0,Q0)
C--
C--   SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
C--  &                   PHI0,FMASK,TLAND,TSEA,SWAV,SSRD,SLRD,
C--  &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
C--  &                   TSFC,TSKIN,U0,V0,T0,Q0)
C--
C--   Purpose: Compute surface fluxes of momentum, energy and moisture,
C--            and define surface skin temperature from energy balance
C--   Input:   PSA    = norm. surface pressure [p/p0]   (2-dim)
C--            UA     = u-wind                          (3-dim)
C--            VA     = v-wind                          (3-dim)
C--            TA     = temperature                     (3-dim)
C--            QA     = specific humidity [g/kg]        (3-dim)
C--            RH     = relative humidity [0-1]         (3-dim)
C--            PHI    = geopotential                    (3-dim)
C--            PHI0   = surface geopotential            (2-dim)
C--            FMASK  = fractional land-sea mask        (2-dim)
C--            TLAND  = land-surface temperature        (2-dim)
C--            TSEA   =  sea-surface temperature        (2-dim)
C--            SWAV   = soil wetness availability [0-1] (2-dim)
C--            SSRD   = sfc sw radiation (downw. flux)  (2-dim)
C--            SLRD   = sfc lw radiation (downw. flux)  (2-dim)
C--   Output:  USTR   = u stress                        (2-dim)
C--            VSTR   = v stress                        (2-dim)
C--            SHF    = sensible heat flux              (2-dim)
C--            EVAP   = evaporation [g/(m^2 s)]         (2-dim)
C--            SLRU   = sfc lw radiation (upward flux)  (2-dim)
C--            HFLUXN = net heat flux into land/sea     (2-dim)           
C--            TSFC   = surface temperature (clim.)     (2-dim)
C--            TSKIN  = skin surface temperature        (2-dim)
C--            U0     = near-surface u-wind             (2-dim)
C--            V0     = near-surface v-wind             (2-dim)
C--            T0     = near-surface air temperature    (2-dim)
C--            Q0     = near-surface sp. humidity [g/kg](2-dim)
C--
C--IO h atparam.h, atparam1.h
C--IO h com_physcon.h, com_sflcon.h, com_radcon.h
C--IO h planetparam.h, com_planet.h
C--IO sx 288. - temperature in degrees Kelvin - TTROP
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"
      include "planetparam.h"

C     Physical constants + functions of sigma and latitude

      include "com_planet.h"
      include "com_physcon.h"

C     Surface flux constants

      include "com_sflcon.h"      

      include "com_radcon.h"

      REAL PSA(NGP), UA(NGP,NLEV), VA(NGP,NLEV), TA(NGP,NLEV),
     &     QA(NGP,NLEV), RH(NGP,NLEV), PHI(NGP,NLEV),
     &     PHI0(NGP), FMASK(NGP), TLAND(NGP), TSEA(NGP), SWAV(NGP),
     &     SSRD(NGP), SLRD(NGP)

      REAL USTR(NGP,3), VSTR(NGP,3), SHF(NGP,3), EVAP(NGP,3),
     &     SLRU(NGP,3), HFLUXN(NGP,2), TSFC(NGP), TSKIN(NGP),
     &     U0(NGP), V0(NGP), T0(NGP), Q0(NGP)

      REAL T1(NGP,2), T2(NGP,2), Q1(NGP,2), QSAT0(NGP,2), 
     &     DENVVS(NGP,0:2), DSLR(NGP), DTSKIN(NGP), CLAMB(NGP)

      LOGICAL LSCASYM, LSCDRAG, LSKINEB

      LOGICAL LFLUXLAND

      SAVE T1, Q1, DENVVS

      LSCASYM = .true.   ! true : use an asymmetric stability coefficient
      LSCDRAG = .true.   ! true : use stability coef. to compute drag over sea
      LSKINEB = .true.   ! true : redefine skin temp. from energy balance
  
c      CLAMBDA = 7.       ! Heat conductivity in skin layer
c      CLAMBSN = 7.       ! Heat conductivity for snow cover = 1

      ESBC  = EMISFC*SBC
      ESBC4 = 4.*ESBC

      GHUM0 = 1.-FHUM0
 
      DLAMBDA = CLAMBSN-CLAMBDA

C--   1. Extrapolation of wind, temp, hum. and density to the surface

C     1.1 Wind components
   
      DO J=1,NGP
        U0(J) = FWIND0*UA(J,NLEV)
        V0(J) = FWIND0*VA(J,NLEV)
      ENDDO

C     1.2 Temperature

      GTEMP0 = 1.-FTEMP0
      RCP = 1./CP
      RDPHI0 =-1./(RD*TTROP*SIGL(NLEV))
      NL1=NLEV-1

      DO J=1,NGP
c       Temperature difference between lowest level and sfc
        DT1 = WVI(NLEV,2)*(TA(J,NLEV)-TA(J,NL1))
c       Extrapolated temperature using actual lapse rate (1:land, 2:sea)
        T1(J,1) = TA(J,NLEV)+DT1
        T1(J,2) = T1(J,1)+PHI0(J)*DT1*RDPHI0
c       Extrapolated temperature using dry-adiab. lapse rate (1:land, 2:sea)
        T2(J,2) = TA(J,NLEV)+RCP*PHI(J,NLEV)
        T2(J,1) = T2(J,2)-RCP*PHI0(J)
      ENDDO

      DO J=1,NGP
        IF (TA(J,NLEV).GT.TA(J,NL1)) THEN
C         Use extrapolated temp. if dT/dz < 0
          T1(J,1) = FTEMP0*T1(J,1)+GTEMP0*T2(J,1)
          T1(J,2) = FTEMP0*T1(J,2)+GTEMP0*T2(J,2)
        ELSE
C         Use temp. at lowest level if dT/dz > 0
          T1(J,1) = TA(J,NLEV)
          T1(J,2) = TA(J,NLEV)
        ENDIF
        T0(J) = T1(J,2)+FMASK(J)*(T1(J,1)-T1(J,2))
      ENDDO

C     1.3 Spec. humidity

C      GHUM0 = 1.-FHUM0

c      CALL SHTORH (-1,NGP,T0,PSA,1.,Q0,RH(1,NLEV),QSAT0)

c      DO J=1,NGP
c        Q0(J)=FHUM0*Q0(J)+GHUM0*QA(J,NLEV)
c      ENDDO

C     1.3 Density * wind speed (including gustiness factor)

      PRD = P0/RD
      VG2 = VGUST*VGUST

      DO J=1,NGP
        DENVVS(J,0)=(PRD*PSA(J)/T0(J))*
     &              SQRT(U0(J)*U0(J)+V0(J)*V0(J)+VG2)
      ENDDO

C     2. Compute land-sfc. fluxes using prescribed skin temperature

C     2.1 Define effective skin temperature to compensate for
C         non-linearity of heat/moisture fluxes during the daily cycle

C     2.2 Stability correction = f[pot.temp.(sfc)-pot.temp.(air)]  

      RDTH  = FSTAB/DTHETA
      ASTAB = 1.
      IF (LSCASYM) ASTAB = 0.5   ! to get smaller dS/dT in stable conditions

      DO J=1,NGP
        IF (TSEA(J).GT.T2(J,2)) THEN
           DTHS=MIN(DTHETA,TSEA(J)-T2(J,2))
        ELSE
           DTHS=MAX(-DTHETA,ASTAB*(TSEA(J)-T2(J,2)))
        ENDIF
        DENVVS(J,2)=DENVVS(J,0)*(1.+DTHS*RDTH)
        TSKIN(J)=TSEA(J)
      ENDDO


      IF (FHUM0.GT.0.) THEN

        CALL SHTORH(-1,NGP,T1(1,2),PSA,1.,Q1(1,2),RH(1,NLEV),QSAT0(1,2))

        DO J=1,NGP
          Q1(J,2) = FHUM0*Q1(J,2)+GHUM0*QA(J,NLEV)
        ENDDO

      ELSE

        DO J=1,NGP
          Q1(J,2) = QA(J,NLEV)
        ENDDO

      ENDIF

C     4.2 Wind stress

      KS=2
      IF (LSCDRAG) KS = 2

      DO J=1,NGP
        CDSDV     =  CDS*DENVVS(J,KS)
        USTR(J,2) = -CDSDV*UA(J,NLEV)
        VSTR(J,2) = -CDSDV*VA(J,NLEV)
      ENDDO


C     Start of sea-sfc. heat fluxes computation

C     4.3 Sensible heat flux 

      KS=2
      CHSCP = CHS*CP

      DO J=1,NGP
        SHF(J,2) = CHSCP*DENVVS(J,KS)*(TSEA(J)-T1(J,2))
      ENDDO

C     4.4 Evaporation

      CALL SHTORH (0,NGP,TSEA,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,2))

      DO J=1,NGP
        EVAP(J,2) = CHS*DENVVS(J,KS)*(QSAT0(J,2)-Q1(J,2))
      ENDDO

C     4.5 Emission of lw radiation from the surface
C         and net heat fluxes into sea surface

      DO J=1,NGP
        SLRU(J,2)   = ESBC*TSEA(J)**4
        HFLUXN(J,2) = SSRD(J)*(1.-ALB_S(J))+SLRD(J)-
     &                (SLRU(J,2)+SHF(J,2)+ALHC*EVAP(J,2))
      ENDDO

C     End of sea-sfc. heat fluxes computation

C--   3. Weighted average of surface fluxes and temperatures 
C--      according to land-sea mask

        DO J=1,NGP
          USTR(J,3) = USTR(J,2)
          VSTR(J,3) = VSTR(J,2)
           SHF(J,3) =  SHF(J,2)
          EVAP(J,3) = EVAP(J,2)
          SLRU(J,3) = SLRU(J,2)
          USTR(J,1) = USTR(J,2)
          VSTR(J,1) = VSTR(J,2)
           SHF(J,1) =  SHF(J,2)
          EVAP(J,1) = EVAP(J,2)
          SLRU(J,1) = SLRU(J,2)
        ENDDO

        DO J=1,NGP
          TSFC(J)  = TSEA(J)
          TSKIN(J) = TSEA(J)
          T0(J)    = T1(J,2)
          Q0(J)    = Q1(J,2)
        ENDDO

      RETURN
      END

      SUBROUTINE SFLSET (PHI0)
C--
C--   SUBROUTINE SFLSET (PHI0)
C--
C--   Purpose: compute orographic factor for land surface drag
C--   Input:   PHI0   = surface geopotential            (2-dim)
C--            Initialized common blocks: SFLFIX
C--IO h atparam.h, atparam1.h, com_physcon.h, com_sflcon.h

C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"

C     Physical constants + functions of sigma and latitude
      include "com_physcon.h"

C     Surface flux constants
      include "com_sflcon.h"

      REAL PHI0(NLON,NLAT)

      RHDRAG = 1./(GG*HDRAG)

      DO I=1,NLON
        DO J=1,NLAT
          FOROG(I,J)=1.+FHDRAG*(1.-EXP(-MAX(PHI0(I,J),0.)*RHDRAG))
        ENDDO
      ENDDO

C--
      RETURN
      END
