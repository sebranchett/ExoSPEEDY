
      SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
     &                   PHI0,FMASK,TLAND,TSEA,SWAV,SSRD,SLRD,
     &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
     &                   TSFC,TSKIN,U0,V0,T0,Q0,LFLUXLAND)
C--
C--   SUBROUTINE SUFLUX (PSA,UA,VA,TA,QA,RH,PHI,
C--  &                   PHI0,FMASK,TLAND,TSEA,SWAV,SSRD,SLRD,
C--  &                   USTR,VSTR,SHF,EVAP,SLRU,HFLUXN,
C--  &                   TSFC,TSKIN,U0,V0,T0,Q0,LFLUXLAND)
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
C--            LFLUXLAND   = Logical related ti flux-correction
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
C     Resolution parameters

      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude

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

      IF ( LFLUXLAND )  THEN

C--   1. Extrapolation of wind, temp, hum. and density to the surface

C     1.1 Wind components
   
      DO J=1,NGP
        U0(J) = FWIND0*UA(J,NLEV)
        V0(J) = FWIND0*VA(J,NLEV)
      ENDDO

C     1.2 Temperature

      GTEMP0 = 1.-FTEMP0
      RCP = 1./CP
      RDPHI0 =-1./(RD*288.*SIGL(NLEV))
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

      DO JLAT=1,NLAT
	J0=NLON*(JLAT-1)
        SQCLAT=SQRT(CLAT(JLAT))
        DO J=J0+1,J0+NLON
          TSKIN(J)=TLAND(J)+CTDAY*SQCLAT*SSRD(J)*(1.-ALB_L(J))*PSA(J)
        ENDDO
      ENDDO

C     2.2 Stability correction = f[pot.temp.(sfc)-pot.temp.(air)]  

      RDTH  = FSTAB/DTHETA
      ASTAB = 1.
      IF (LSCASYM) ASTAB = 0.5   ! to get smaller dS/dT in stable conditions

      DO J=1,NGP

C       potential temp. difference (land+sea average)
cfk        DTH0 = TSEA(J)-T2(J,2)
cfk        DTH0 = DTH0+FMASK(J)*((TSKIN(J)-T2(J,1))-DTH0)

cfk        IF (DTH0.GT.0.0) THEN
cfk           DTHL=MIN(DTHETA,DTH0)
cfk        ELSE
cfk           DTHL=MAX(-DTHETA,ASTAB*DTH0)
cfk        ENDIF

cfk        DENVVS(J,1)=DENVVS(J,0)*(1.+DTHL*RDTH)

        IF (TSKIN(J).GT.T2(J,1)) THEN
           DTHL=MIN(DTHETA,TSKIN(J)-T2(J,1))
        ELSE
           DTHL=MAX(-DTHETA,ASTAB*(TSKIN(J)-T2(J,1)))
        ENDIF
        DENVVS(J,1)=DENVVS(J,0)*(1.+DTHL*RDTH)

      ENDDO

C     2.3 Wind stress 

      DO J=1,NGP
        CDLDV     =  CDL*DENVVS(J,0)*FOROG(J)
        USTR(J,1) = -CDLDV*UA(J,NLEV)
        VSTR(J,1) = -CDLDV*VA(J,NLEV)
      ENDDO

C     2.4 Sensible heat flux 

      CHLCP = CHL*CP

      DO J=1,NGP
        SHF(J,1) = CHLCP*DENVVS(J,1)*(TSKIN(J)-T1(J,1))
      ENDDO

C     2.5 Evaporation

      IF (FHUM0.GT.0.) THEN

        CALL SHTORH(-1,NGP,T1(1,1),PSA,1.,Q1(1,1),RH(1,NLEV),QSAT0(1,1))

        DO J=1,NGP
          Q1(J,1) = FHUM0*Q1(J,1)+GHUM0*QA(J,NLEV)
        ENDDO

      ELSE

        DO J=1,NGP
          Q1(J,1) = QA(J,NLEV)
        ENDDO

      ENDIF

      CALL SHTORH (0,NGP,TSKIN,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,1))

      DO J=1,NGP
C       EVAP(J,1) = CHL*DENVVS(J,1)*SWAV(J)*MAX(0.,QSAT0(J,1)-Q1(J,1))
        EVAP(J,1) = CHL*DENVVS(J,1)*MAX(0.,SWAV(J)*QSAT0(J,1)-Q1(J,1))
      ENDDO

C--   3. Compute land-surface energy balance;
C--      adjust skin temperature and heat fluxes

C     3.1. Emission of lw radiation from the surface
C          and net heat fluxes into land surface

      DO J=1,NGP
        TSK3        = TSKIN(J)**3
        DSLR(J)     = ESBC4*TSK3
        SLRU(J,1)   = ESBC *TSK3*TSKIN(J)
        HFLUXN(J,1) = SSRD(J)*(1.-ALB_L(J))+SLRD(J)-
     &                (SLRU(J,1)+SHF(J,1)+ALHC*EVAP(J,1))
      ENDDO

C     3.2 Re-definition of skin temperature from energy balance

      IF ( LSKINEB ) THEN

C       Compute net heat flux including flux into ground
        DO J=1,NGP
          CLAMB(J)    = CLAMBDA+SNOWC(J)*DLAMBDA
          HFLUXN(J,1) = HFLUXN(J,1)-CLAMB(J)*(TSKIN(J)-TLAND(J))
          DTSKIN(J)   = TSKIN(J)+1.
        ENDDO

C       Compute d(Evap) for a 1-degree increment of Tskin

        CALL SHTORH (0,NGP,DTSKIN,PSA,1.,QDUMMY,RDUMMY,QSAT0(1,2))

        DO J=1,NGP
          IF (EVAP(J,1).GT.0) THEN
             QSAT0(J,2) = SWAV(J)*(QSAT0(J,2)-QSAT0(J,1))
          ELSE
             QSAT0(J,2) = 0.
          ENDIF
        ENDDO

C       Redefine skin temperature to balance the heat budget 
        DO J=1,NGP
          DHFDT     = CLAMB(J)+DSLR(J)+
     &                CHL*DENVVS(J,1)*(CP+ALHC*QSAT0(J,2))
          DTSKIN(J) = HFLUXN(J,1)/DHFDT
          TSKIN(J)  = TSKIN(J)+DTSKIN(J)
        ENDDO

C       Add linear corrections to heat fluxes
        DO J=1,NGP
          SHF(J,1)    = SHF(J,1) +CHLCP*DENVVS(J,1)*DTSKIN(J)
          EVAP(J,1)   = EVAP(J,1)+CHL*DENVVS(J,1)*QSAT0(J,2)*DTSKIN(J)
          SLRU(J,1)   = SLRU(J,1)+DSLR(J)*DTSKIN(J)
          HFLUXN(J,1) = CLAMB(J)*(TSKIN(J)-TLAND(J))
        ENDDO

      ENDIF

c      ENDIF

C--   4. Compute sea surface fluxes:
C--      Note: stability terms and wind stress are NOT re-defined
C--            if LFLUXLAND = .false.

C     4.1 Correct near-sfc. air temperature over coastal sea points
C         and compute near-sfc. humidity
 

cfk      DO J=1,NGP
cfk         if (FMASK(j).gt.0.) then
cfk            dtsea  = TSEA(J) -T1(J,2)
cfk            dtland = TSKIN(J)-T1(J,1)
cfk            if (dtsea.gt.0.0.and.dtland.lt.0.0) then
cfk               dtsea   = dtsea*(1.-FMASK(j)**2)
cfk               T1(J,2) = TSEA(J)-dtsea
cfk            endif
cfk         endif
cfk      ENDDO


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

cfk      KS = 0
      KS=2
cfk      IF (LSCDRAG) KS = 1
      IF (LSCDRAG) KS = 2

      DO J=1,NGP
        CDSDV     =  CDS*DENVVS(J,KS)
        USTR(J,2) = -CDSDV*UA(J,NLEV)
        VSTR(J,2) = -CDSDV*VA(J,NLEV)
      ENDDO

C     End of 'land-mode' computation
      ENDIF

C     Start of sea-sfc. heat fluxes computation

C     4.3 Sensible heat flux 

cfk      KS = 1
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

      IF ( LFLUXLAND )  THEN

        DO J=1,NGP
          USTR(J,3) = USTR(J,2)+FMASK(J)*(USTR(J,1)-USTR(J,2))
          VSTR(J,3) = VSTR(J,2)+FMASK(J)*(VSTR(J,1)-VSTR(J,2))
           SHF(J,3) =  SHF(J,2)+FMASK(J)*( SHF(J,1)- SHF(J,2))
          EVAP(J,3) = EVAP(J,2)+FMASK(J)*(EVAP(J,1)-EVAP(J,2))
          SLRU(J,3) = SLRU(J,2)+FMASK(J)*(SLRU(J,1)-SLRU(J,2))
        ENDDO

        DO J=1,NGP
          TSFC(J)  = TSEA(J)+FMASK(J)*(TLAND(J)-TSEA(J))
          TSKIN(J) = TSEA(J)+FMASK(J)*(TSKIN(J)-TSEA(J))
          T0(J)    = T1(J,2)+FMASK(J)*(T1(J,1)- T1(J,2))
          Q0(J)    = Q1(J,2)+FMASK(J)*(Q1(J,1)- Q1(J,2))
        ENDDO

      ENDIF

      RETURN
      END

      SUBROUTINE SFLSET (PHI0)
C--
C--   SUBROUTINE SFLSET (PHI0)
C--
C--   Purpose: compute orographic factor for land surface drag
C--   Input:   PHI0   = surface geopotential            (2-dim)
C--            Initialized common blocks: SFLFIX

C     Resolution parameters
C
      include "atparam.h"
      include "atparam1.h"

      PARAMETER ( NLON=IX, NLAT=IL, NLEV=KX, NGP=NLON*NLAT )

C     Physical constants + functions of sigma and latitude
      include "com_physcon.h"

C     Surface flux constants
      include "com_sflcon.h"

      REAL PHI0(NGP)

      RHDRAG = 1./(GG*HDRAG)

      DO J=1,NGP
        FOROG(J)=1.+FHDRAG*(1.-EXP(-MAX(PHI0(J),0.)*RHDRAG))
      ENDDO

C--
      RETURN
      END
