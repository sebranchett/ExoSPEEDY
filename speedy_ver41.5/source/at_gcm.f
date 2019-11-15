      PROGRAM AGCM_MAIN
C--
C--   Program : AGCM_MAIN
C--

      CHARACTER (len=3) cexp        ! experiment identifier

C--   1. Initialization
C--      ndays = no. of integration days, set by AGCM_INIT
 
      cexp='exp'
      CALL AGCM_INIT (cexp,0,0,0,
     &                 ndays)

      print *, 'Integration length in days: ', ndays

C--   2. Do loop over total no. of integration days

      DO jday = 1,ndays

C--      2.2 Run atmospheric model for 1 day

         CALL AGCM_1DAY (jday)

C--      2.1 Exchange data with coupler

         CALL AGCM_TO_COUPLER (jday)

         CALL COUPLER_TO_AGCM (jday)

C--      2.3 Write restart file at the end of selected months
C--          and at the end of the integration 

         CALL RESTART (jday)

      ENDDO

      STOP
      END


      SUBROUTINE AGCM_1DAY (jday)
C--
C--   SUBROUTINE AGCM_1DAY (jday)
C--
C--   Perform atm. model integration for 1 day, 
C--   post-proc. and I/O at selected times 
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_tsteps.h"
      include "com_date.h"


      if (iday.eq.1) PRINT *, ' Start of year/month = ', iyear, imonth

      istep = 1 + (jday-1)*nsteps

C--   1. Set forcing terms according to date

      CALL FORDATE (1)

C--   2. Set daily-average flux arrays to zero

      CALL DMFLUX (0)

C--   3. Integrate the atmospheric model for 1 day

      CALL STLOOP (istep)

C--   4. Write daily-mean output

      CALL DMOUT (idout)

C--   5. Compute new date

      CALL NEWDATE (1)

C--   6. Write time-mean output files and restart file
C--      at the end of selected months

      if (iday.eq.1) then

C         Write monthly-mean output for previous month
          if (nstout.lt.0) CALL TMOUT (1)

C         Open new output files at the beginning of each year
          if (imonth.eq.1.and.jday.lt.ndaytot) CALL SETGRD (1,cexp)

      endif

      RETURN
      END
