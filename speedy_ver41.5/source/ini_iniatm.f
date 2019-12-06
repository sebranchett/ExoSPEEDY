      SUBROUTINE INI_ATM (cexp)
C--
C--   SUBROUTINE INI_ATM (cexp)
C--
C--   Purpose : Call initialization routines for all model common blocks 
C--
C--IO h atparam.h, atparam1.h
C--IO h com_tsteps.h, com_date.h, par_tmean.h, com_dyncon1.h
C--IO s 1440 = minutes in a day?
C--IO w use SETCTL to write 'attm' to unit 12
C--IO w use SETCTL to write 'atva' to unit 14
C--IO w use SETCTL to write 'atdf' to unit 16
C--IO w use SETCTL_D to write 'daytm' to unit 18
C--IO w use SETGRD to write year to units 11, 13 and 15
C--IO w use TMOUT to write time-means and variances to units 11, 13 and 15
C--IO w use DMOUT to write daily-means to unit 17

      include "atparam.h"
      include "atparam1.h"

      include "com_tsteps.h"
      include "com_date.h"

      include "par_tmean.h"

      include "com_dyncon1.h"

      CHARACTER*3 cexp        ! experiment identifier

      REAL PPL(KX)            ! post-processing levels (hPa/1000)

      iitest=1

C--   1. Initialize ffts

      if (iitest.eq.1) print*, 'calling INIFFT'
      CALL INIFFT

C--   2. Initialize dynamical constants and operators
 
      if (iitest.eq.1) print*, 'calling INDYNS'
      CALL INDYNS

C--   3. Set post-processing levels

      do k=1,kx
        PPL(k)=PRLEV(FSG(k))
      enddo

C--   4. Initialize constants for physical parametrization

      if (iitest.eq.1) print*, 'calling INPHYS'
      CALL INPHYS (HSG,PPL,RADANG)

C--   5. Initialize forcing fields (boundary cond. + random forcing)

      if (iitest.eq.1) print*, 'calling INBCON'
      CALL INBCON (grav,RADANG)

      if (iitest.eq.1) print*, 'calling INIRDF'
      CALL INIRDF (indrdf)

C--   6. Initialize model variables

      if (iitest.eq.1) print*, 'calling INVARS'
      CALL INVARS

C--   7. Initialize time-mean arrays for surface fluxes and output fields

      if (iitest.eq.1) print*, 'calling DMFLUX'
      CALL DMFLUX (0)

      if (iitest.eq.1) print*, 'calling DMOUT'
      CALL DMOUT (0)

      if (iitest.eq.1) print*, 'calling TMOUT'
      CALL TMOUT (0)
 
C--   8. Set up the time-mean and daily-mean output (grads format)

C     8.1 Control files for time-means

      if (nstout.le.0) then
         ntm  = nmonts
         ndtm = -1
      else
         ntm  = ndaytot*nsteps/nstout
         ndtm = 1440*nstout/nsteps
      endif

      if (iitest.eq.1) print*, 'calling SETCTL'

      IS3D=1
      CALL SETCTL (12,IX,IL,KX,NTM,NDTM,IS3D,NS3D1,NS2D_1,NS2D_2,
     *             RADANG,PPL,'attm',cexp,IYEAR0,IMONT0)
      IS3D=IS3D+NS3D1
      CALL SETCTL (14,IX,IL,KX,NTM,NDTM,IS3D,NS3D2,0,0,
     *             RADANG,PPL,'atva',cexp,IYEAR0,IMONT0)

      IS3D=IS3D+NS3D2
      CALL SETCTL (16,IX,IL,KX,NTM,NDTM,IS3D,NS3D3,0,0,
     *             RADANG,PPL,'atdf',cexp,IYEAR0,IMONT0)

C     8.2 Control files for daily means

      ndm =  ndaytot
      nddm = 1

      if (iitest.eq.1) print*, 'calling SETCTL_D'

      if (IDOUT.eq.1) then
         CALL SETCTL_D (18,IX,IL,KX,NDM,NDDM,3,1,
     *                  RADANG,PPL,'daytm',cexp,IYEAR0,IMONT0)
      else if (IDOUT.eq.2) then
         CALL SETCTL_D (18,IX,IL,KX,NDM,NDDM,NS2D_D1,1,
     *                  RADANG,PPL,'daytm',cexp,IYEAR0,IMONT0)
      else if (IDOUT.ge.3) then
         CALL SETCTL_D (18,IX,IL,KX,NDM,NDDM,NS2D_D1,NS2D_D2,
     *                  RADANG,PPL,'daytm',cexp,IYEAR0,IMONT0)
      endif

C     8.3 Output files for grid-point fields

      CALL SETGRD (0,cexp)

C--
      RETURN
      END
  

      FUNCTION PRLEV (SIGLEV)
C--									
C--   FUNCTION PRLEV (SIGLEV)
C--   Purpose : select the closest standard pressure level for post-proc.
C--   Input :   SIGLEV = sigma level
C--IO s PLEV = standard pressure level for post-proc.

      REAL PLEV(14)

      DATA PLEV/ 0.925, 0.850, 0.775, 0.700, 0.600, 0.500, 0.400,  
     *           0.300, 0.250, 0.200, 0.150, 0.100, 0.050, 0.030/

      PRLEV = 1.
      DIF = 1.-SIGLEV

      DO K=1,14
        ADIF = ABS(PLEV(K)-SIGLEV)
        IF (ADIF.LE.DIF) THEN
          DIF = ADIF
          PRLEV = PLEV(K)
        ENDIF
      ENDDO

      RETURN
      END

