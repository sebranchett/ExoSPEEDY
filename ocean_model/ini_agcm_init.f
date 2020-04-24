
      SUBROUTINE AGCM_INIT (cexp,inidate,ntimes,irstart,
     &                      ndays)
C--
C--   SUBROUTINE AGCM_INIT (cexp,inidate,ntimes,irstart,
C--  &                      ndays)
C--
C--   Purpose: Initialization of atmos. model and coupling interface 
C--
c      include "atparam.h"
c      include "atparam1.h"

c      parameter ( NGP = IX*IL )
C--IO h com_tsteps.h, com_date.h, com_lflags.h, com_cpl_flags.h
C--IO h planetparam.h, com_planet.h
C--IO sx concept of months and days

      include "planetparam.h"
      include "com_planet.h"
      include "com_tsteps.h"
      include "com_date.h"

      include "com_lflags.h"
      include "com_cpl_flags.h"


C     Input (reset by input/include files if inidate = 0):
      CHARACTER*3 cexp        ! experiment identifier
      INTEGER     inidate     ! initial date YYYYMM
      INTEGER     ntimes      ! integr. length in months (< 0) or days (> 0)
      INTEGER     irstart     ! restart flag: 0 = no, > 0 = yes

C     Output:
      INTEGER     ndays       ! total no. of integration days

      PRINT *, ' Hallo from Speedy_AGCM'

C--   1. Set run initial time, duration, time-stepping and coupling options

C--IO h cls_instep.h
C--IO h cls_inplanet.h
C--IO r read ISTART, restart flag: 0 = no, > 0 = yes, from unit 2
C--IO r read CEXP, experiment identifier, from unit 2
C--IO sx 12 months in a year
      if (inidate.le.0) then
         read (2,*) istart
         read (2,'(a3)') cexp
      else
         istart = irstart
      endif

      if (istart.ne.0) istart = 1

      include "cls_inplanet.h"
      include "cls_instep.h"

      if (inidate.gt.0) then

         iyear0 = inidate/100
         imont0 = mod(inidate,100)

         isst0 = (iyear0-issty0)*months+imont0

         if (ntimes.lt.0) then
            nmonts = -ntimes
            ndaysl = 0
         else
            nmonts = 0
            ndaysl = ntimes
         endif

      endif

      CALL NEWDATE (0)

      ndays = ndaytot

C     check consistency of coupling and prescribed SST anomaly flags
      if (icsea.ge.4) isstan = 1

C--   2. Initialization of atmospheric model constants and variables 

      CALL INI_ATM (cexp)

C--   3. Initialization of coupled modules (land, sea, ice)
 
      CALL INI_COUPLER (istart)

C--   4. Set up the forcing fields for the first time step

      CALL FORDATE (0)

C--   5. Do the initial (2nd-order) time step, initialize the semi-impl. scheme

      CALL STEPONE

C--
      RETURN
      END


      subroutine NEWDATE (imode)
C--
C--   subroutine NEWDATE (imode)
C--   Purpose:   initilialize and update date variables 
C--   Input :    imode = 0 for initialization, > 0 for update  
C--IO h com_date.h, com_tsteps.h
C--IO h planetparam.h, com_planet.h
C--IO sx 365 days in a year
C--IO sx 12 months in a year
C--IO sx days in each month: 31,28,31,30,31,30,31,31,30,31,30,31

      parameter ( NCAL=360 )

      include "planetparam.h"
      include "com_planet.h"
      include "com_date.h"
      include "com_tsteps.h"

      if (imode.le.0) then

C        Calendar

         ndaycal(:,1) = 30
 
         ndaycal(1,2) = 0
         if (months.gt.1) then
            do jm=2,months
               ndaycal(jm,2) = ndaycal(jm-1,1)+ndaycal(jm-1,2)
            enddo
         endif

C        Total no. of integration days

         ndaytot = ndaysl
         im = imont0
      
         do jm=1,nmonts
            ndaytot = ndaytot+ndaycal(im,1)
            im = im+1
            if (im.eq.(months+1)) im=1
         enddo

C        Initial date

         iyear  = iyear0
         imonth = imont0
         iday   = 1

      else

C        Set new date

         iday = iday+1

         if (iday.gt.ndaycal(imonth,1)) then
            iday   = 1
            imonth = imonth+1
         endif

         if (imonth.gt.months) then
            imonth = 1
            iyear  = iyear+1
         endif

      endif

C     Additional variables to define forcing terms and boundary cond.

      if (iseasc.ge.1) then

         imont1 = imonth
         tmonth = (iday-0.5)/float(ndaycal(imonth,1))
         tyear  = (ndaycal(imonth,2)+iday-0.5)/float(ncal)

      else

         imont1 = imont0
         tmonth = 0.5
         tyear  = (ndaycal(imont1,2)+0.5*ndaycal(imont1,2))/
     &            float(ncal)

      endif

      return
      end



