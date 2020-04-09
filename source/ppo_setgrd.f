
      SUBROUTINE SETGRD (ind,norun)
C--
C--   SUBROUTINE SETGRD (ind)
C--   Purpose : open and close output files (.grd)
C--
C--   Input : ind = 0 for initialization, 1 otherwise
C--         : norun = run identifier
C--
C--IO h com_date.h, com_tsteps.h
C--IO h planetparam.h, com_lflags.h
C--IO sx .grd file unit numbers 11(attm), 13(atva), 15(atdf) and 17(daytm)

      include "planetparam.h"
      include "com_date.h"
      include "com_tsteps.h"
      include "com_lflags.h"

      character*3  norun
      character*16 ofile11, ofile13, ofile15
      character*17 ofile17
      character*16 ofile91, ofile93, ofile95
      character*17 ofile97

      save ofile91, ofile93, ofile95, ofile97

      save ofile11, ofile13, ofile15, ofile17

      if (ind.eq.0) then

         ofile11='attmNNN_YYYY.grd'
         ofile13='atvaNNN_YYYY.grd'
         ofile15='atdfNNN_YYYY.grd'

         ofile11(5:7)=norun
         ofile13(5:7)=norun
         ofile15(5:7)=norun

         if (IDOUT .gt. 0) then
            ofile17='daytmNNN_YYYY.grd'
            ofile17(6:8)=norun
         endif

      endif

      write (ofile11(9:12),'(i4)') iyear
      write (ofile13(9:12),'(i4)') iyear
      write (ofile15(9:12),'(i4)') iyear

      if (IDOUT.gt.0) write (ofile17(10:13),'(i4)') iyear

      if (LTXTO) then
         ofile91 = ofile11
         ofile93 = ofile13
         ofile95 = ofile15
         ofile97 = ofile17

         write (ofile91(14:16),'(a)') 'txt'
         write (ofile93(14:16),'(a)') 'txt'
         write (ofile95(14:16),'(a)') 'txt'
         if (IDOUT.gt.0) write (ofile97(15:17),'(a)') iyear
      endif

      if (ind.ne.0) then

         close( unit=11 )
         close( unit=13 )
         close( unit=15 )

         if (IDOUT.gt.0) close( unit=17 )

         if (LTXTO) then
            close( unit=91 )
            close( unit=93 )
            close( unit=95 )

            if (IDOUT.gt.0) close( unit=97 )
         endif

      endif
      
      open ( unit=11, file=ofile11, 
     &       status='new', form='unformatted', access='sequential' )
      open ( unit=13, file=ofile13, 
     &       status='new', form='unformatted', access='sequential' )
      open ( unit=15, file=ofile15, 
     &       status='new', form='unformatted', access='sequential' )

      if (IDOUT.gt.0) then
         open ( unit=17, file=ofile17,
     &          status='new', form='unformatted', access='sequential' )
      endif

      if (LTXTO) then
         open ( unit=91, file=ofile91, form='formatted' )
         open ( unit=93, file=ofile93, form='formatted' )
         open ( unit=95, file=ofile95, form='formatted' )

         if (IDOUT.gt.0) then
            open ( unit=97, file=ofile97, form='formatted' )
         endif
      endif

      return
      end


