      SUBROUTINE FORDATE (imode)
C--
C--   SUBROUTINE FORDATE (imode)
C--   
C--   Purpose :	Compute forcing fields for the current date
C--             and correction terms for horiz. diffusion
C--
C--   Input : imode : 0 = initialization step, 1 = daily update

C--IO h atparam.h, atparam1.h, com_tsteps.h, com_date.h
C--IO h com_dyncon0.h, com_physcon.h, com_radcon.h, com_hdifcon.h
C--IO h com_surfcon.h, com_cli_sea.h, com_cli_land.h
C--IO h com_var_sea.h, com_var_land.h, com_lflags.h
C--IO h planetparam.h
C--IO sx RYRCO2 = 1950, CO2 absorptivity
C--IO sx DELCO2 = 0.005, CO2 absorptivity, rate of change per year
C--IO w write PHIS0, surface geopotential, to unit 19
C--IO w write CORH, temperature correction term horiz. diffusion, to unit 19
C--IO w write CORH, Humidity correction term horiz. diffusion, to unit 19

      include "atparam.h"
      include "atparam1.h"
      include "planetparam.h"
C
      include "com_tsteps.h"
      include "com_date.h"

      include "com_dyncon0.h"
      include "com_physcon.h"

      include "com_radcon.h"
      include "com_hdifcon.h"

      include "com_surfcon.h"

      include "com_cli_sea.h"
      include "com_cli_land.h"

      include "com_var_sea.h"
      include "com_var_land.h"

      include "com_lflags.h"

      REAL GAMLAT(NLAT),
     &     CORH(NLON,NLAT), TSFC(NLON,NLAT), TREF(NLON,NLAT),
     &     PSFC(NLON,NLAT), QSFC(NLON,NLAT), QREF(NLON,NLAT)

      REAL FLAND(NGP), ALB_0(NGP)
      EQUIVALENCE (FLAND,FMASK_L), (ALB_0,ALB0)

      iitest = 0

C--   Time variables for interpolation are set by NEWDATE

C--   1. Time-independent parts of physical parametrizations

      IF (IMODE.EQ.0) THEN

        CALL RADSET

        CALL SFLSET (PHIS0)

        ABLCO2_ref = ABLCO2

      ENDIF

C--   2. Daily-mean radiative forcing 

C     Incoming solar radiation

cfk      CALL SOL_OZ (TYEAR)

C     Total surface albedo

      do j=1,ngp
         SNOWC(j)  = min(1.,SNOWD_AM(j)/sd2sc)
         ALB_L(j)  = ALB_0(j)+SNOWC(j)*(albsn-ALB_0(j))
         ALB_S(j)  = albsea+SICE_AM(j)*(albice-albsea)
         ALBSFC(j) = ALB_S(j)+FLAND(j)*(ALB_L(j)-ALB_S(j))
      enddo

C     Linear trend of CO2 absorptivity
C     RYRCO2 = reference year
C     DELCO2 = rate of change per year

      IF (LCO2) THEN
         ABLCO2 = ABLCO2_ref*EXP(DELCO2*(IYEAR+TYEAR-RYRCO2))
      ENDIF

C--   3. Temperature correction term for horizontal diffusion

      CALL SETGAM (TYEAR,GAMLAT)

      DO J=1,NLAT
        DO I=1,NLON
          CORH(I,J)=GAMLAT(J)*PHIS0(I,J)
        ENDDO
      ENDDO

      if (iitest.gt.1.and.imode.eq.0) then
         call outest (19,PHIS0)
         call outest (19,CORH)
      endif

      CALL SPEC (CORH,TCORH)

C--   4. Humidity correction term for horizontal diffusion

      ij = 0
      DO J=1,NLAT
        PEXP=1./(RD*GAMLAT(J))
        DO I=1,NLON
          ij = ij+1
c          TSFC(i,j) = FMASK_L(i,j)*STLCL_OB(ij)
c     &               +FMASK_S(i,j)*SSTCL_OB(ij)
          TSFC(i,j) = FMASK_L(i,j)*STL_AM(ij)
     &               +FMASK_S(i,j)*SST_AM(ij)
          TREF(I,J) = TSFC(I,J)+CORH(I,J)
          PSFC(I,J) = (TSFC(I,J)/TREF(I,J))**PEXP
        ENDDO
      ENDDO

      CALL SHTORH (0,NGP,TREF,  1.,-1.,DUMMY,DUMMY,QREF)
      CALL SHTORH (0,NGP,TSFC,PSFC, 1.,DUMMY,DUMMY,QSFC)

      DO J=1,NLAT
        DO I=1,NLON
          CORH(I,J)=REFRH1*(QREF(I,J)-QSFC(I,J))
        ENDDO
      ENDDO

      if (iitest.gt.1.and.imode.eq.0) call outest (19,CORH)

      CALL SPEC (CORH,QCORH)

 900  CONTINUE

      RETURN
      END


      SUBROUTINE SETGAM (TYEAR,GAMLAT)  

C--   Aux. routine GAMLAT : compute reference lapse rate 
C--                         as a function of latitude and date
C--IO h atparam.h, atparam1.h, com_dyncon0.h, com_physcon.h
C--IO sx 1000. conversion from km to meters

      include "atparam.h"
      include "atparam1.h"
C
      include "com_dyncon0.h"
      include "com_physcon.h"

      REAL GAMLAT(NLAT)

      GAMLAT(1) = GAMMA/(1000.*GG)
      DO J=2,NLAT
        GAMLAT(J) = GAMLAT(1)
      ENDDO
C--
      RETURN
      END


      SUBROUTINE OUTEST (iunit,fout)

C--   Aux. routine OUTEST : write one field on a test output file 
C--IO h atparam.h
C--IO w write a matrix to an output file

      include "atparam.h"

      real*4 r4out(ix,il)
      real fout(ix,il)

      do j=1,il
        do i=1,ix
          r4out(i,j)=fout(i,j)
        enddo
      enddo

      write (iunit) r4out

C--
      RETURN
      END

