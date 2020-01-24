      SUBROUTINE INPHYS (HSG,PPL,RLAT)
C--
C--   SUBROUTINE INPHYS (HSG,PPL,RLAT)
C--
C--   Purpose: Initialize common blocks for physical parametrization routines 
C--   Input :  HSG  : sigma at half levels
C--            PPL  : pressure levels for post-processing
C--            RLAT : gaussian-grid latitudes
C--   Initialized common blocks: PHYCON, FSIGLT, FORCON, 
C--                              CNVCON, LSCCON, RADCON, SFLCON, VDICON
C--
C     Resolution parameters
C
C--IO h atparam.h, atparam1.h, com_physcon.h, com_surfcon.h
C--IO h com_cnvcon.h, com_lsccon.h, com_radcon.h, com_sflcon.h
C--IO h com_vdicon.h, cls_inphys.h
C--IO s P0 = 1.e+5
C--IO s GG = 9.81
C--IO s RD = 287.
C--IO s CP = 1004.
C--IO s ALHC = 2501.0
C--IO s        Latent heat is in J/g for consistency with spec.hum. in g/Kg
C--IO s SBC = 5.67e-8
      include "atparam.h"
      include "atparam1.h"
C
C     Physical constants + functions of sigma and latitude
C
      include "com_physcon.h"

C     Surface properties

      include "com_surfcon.h"
C
C     Constants for sub-grid-scale physics
C
      include "com_cnvcon.h"
      include "com_lsccon.h"
      include "com_radcon.h"  
      include "com_sflcon.h"
      include "com_vdicon.h"
C
      REAL HSG(0:NLEV), PPL(NLEV), RLAT(NLAT)  
C
C---  1. Time independent parameters and arrays
C
C     1.1 Physical constants
C
      P0 = 1.e+5
      GG = 9.81
      RD = 287.
      CP = 1004.
C     Latent heat is in J/g for consistency with spec.hum. in g/Kg
      ALHC = 2501.0
cfk#if defined(KNMI)
cfk      ALHS = 2801.0
cfk#endif
      SBC = 5.67e-8
C
C     1.2 Functions of sigma and latitude
C
      SIGH(0)=HSG(0)
C
      DO K=1,NLEV
       SIG(K)  = 0.5*(HSG(K)+HSG(K-1))
       SIGL(K) = LOG(SIG(K))
       SIGH(K) = HSG(K)
       DSIG(K) = HSG(K)-HSG(K-1)
       POUT(K) = PPL(K)
       GRDSIG(K) = GG/(DSIG(K)*P0)
       GRDSCP(K) = GRDSIG(K)/CP
      ENDDO
C
C     Weights for vertical interpolation at half-levels(1,nlev) and surface
C     Note that for phys.par. half-lev(k) is between full-lev k and k+1 
C     Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
C     Fsurf = Ffull(nlev)+WVI(nlev,2)*(Ffull(nlev)-Ffull(nlev-1))
C
      DO K=1,NLEV-1
       WVI(K,1)=1./(SIGL(K+1)-SIGL(K))
       WVI(K,2)=(LOG(SIGH(K))-SIGL(K))*WVI(K,1)
      ENDDO
C
      WVI(NLEV,1)=0.
      WVI(NLEV,2)=(LOG(0.99)-SIGL(NLEV))*WVI(NLEV-1,1)
C
      DO J=1,NLAT
       SLAT(J)=SIN(RLAT(J))
       CLAT(J)=COS(RLAT(J))
      ENDDO

C--   2. Constants for physical parametrization routines:

      include "cls_inphys.h"

C---
      RETURN
      END

