
      SUBROUTINE INI_LAND (istart)
C--
C--   SUBROUTINE INI_LAND (istart)
C-- 
C--   Input : istart = restart flag ( 0 = no, 1 = yes)
C--IO h atparam.h, com_cpl_flags.h, com_cli_land.h, com_var_land.h
C--IO h planetparam.h

      include "atparam.h"
      include "planetparam.h"

      include "com_cpl_flags.h"

      include "com_cli_land.h"
      include "com_var_land.h"

C--   1. Compute climatological fields for initial date

      CALL ATM2LAND (0)

C--   2. Initialize prognostic variables of land model
C--      in case of no restart or no coupling

      if (istart.le.0) then

        stl_lm(:)  = stlcl_ob(:)      ! land sfc. temperature 

      endif

C--   3. Compute additional land variables

      CALL LAND2ATM (0)

      RETURN
      END


      SUBROUTINE ATM2LAND (jday)
C--
C--   SUBROUTINE ATM2LAND (jday)
C-- 
C--IO h atparam.h, com_date.h, com_cpl_flags.h, com_cli_land.h, com_var_land.h
C--IO h com_flx_land.h, com_cplvar_land.h
C--IO h planetparam.h, com_planet.h
      include "atparam.h"
      include "planetparam.h"

      include "com_date.h"
      include "com_cpl_flags.h"

      include "com_planet.h"
      include "com_cli_land.h" 
      include "com_var_land.h"
      include "com_flx_land.h"

      include "com_cplvar_land.h"


C--   1. Interpolate climatological fields to actual date

C     Climatological land sfc. temperature
      call FORIN5 (ngp,imont1,tmonth,stlmnth,stlcl_ob,months)

C     Climatological snow depth
      call FORINT (ngp,imont1,tmonth,snowdmnth,snowdcl_ob,months)

C     Climatological soil water availability
      call FORINT (ngp,imont1,tmonth,soilwmnth,soilwcl_ob,months)

      if (jday.le.0) RETURN

C--   2. Set input variables for mixed-layer/ocean model

      if (icland.gt.0) then

        VLAND_INPUT(:,1) = stl_lm(:)
        VLAND_INPUT(:,2) = hflux_l(:)
        VLAND_INPUT(:,3) = stlcl_ob(:)

      endif

C--   3. Call message-passing routines to send data (if needed)

      RETURN
      END


      SUBROUTINE LAND2ATM (jday)
C--
C--   SUBROUTINE LAND2ATM (jday)
C-- 
C--IO h atparam.h, com_cpl_flags.h, com_var_land.h, com_cplvar_land.h
      include "atparam.h"

      include "com_cpl_flags.h"

      include "com_var_land.h"

      include "com_cplvar_land.h"

      if (jday.gt.0.and.icland.gt.0) then

C--   1. Run ocean mixed layer or 
C--      call message-passing routines to receive data from ocean model

         CALL LAND_MODEL 

C--   2. Get updated variables for mixed-layer/ocean model

         stl_lm(:) = VLAND_OUTPUT(:,1)      ! land sfc. temperature

      endif

C--   3. Compute land-sfc. fields for atm. model

C     3.1 Land sfc. temperature

      if (icland.le.0) then

C        Use observed climatological field
         stl_am(:) = stlcl_ob(:)
         
      else

C        Use land model sfc. temperature
         stl_am(:) = stl_lm(:)

      endif

C     3.2 Snow depth and soil water availability

      snowd_am(:) = snowdcl_ob(:)

      soilw_am(:) = soilwcl_ob(:)


      RETURN
      END


      SUBROUTINE REST_LAND (imode)
C--
C--   SUBROUTINE REST_LAND (imode)
C--
C--   Purpose : read/write land variables from/to a restart file
C--   Input :   IMODE = 0 : read model variables from a restart file
C--                   = 1 : write model variables  to a restart file

C-- 
C--IO h atparam.h, com_cpl_flags.h, com_var_land.h
C--IO r read stl_lm (Land sfc. temperature) from unit (3)
C--IO w write stl_lm, land model variables, to unit (10)
C--IO w write stl_am, atmospheric model fields, to unit (10)
      include "atparam.h"

      include "com_cpl_flags.h"

      include "com_var_land.h"

      if (imode.eq.0) then

         read (3)  stl_lm(:)       ! Land sfc. temperature 

      else

C        write land model variables from coupled runs,
C        otherwise write fields used by atmospheric model

         if (icland.gt.0) then
            write (10) stl_lm(:) 
         else
            write (10) stl_am(:)
         endif

      endif

      RETURN
      END
