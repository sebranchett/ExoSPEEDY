
      SUBROUTINE INI_SEA (istart)
C--
C--   SUBROUTINE INI_SEA (istart)
C-- 
C--   Input : istart = restart flag ( 0 = no, 1 = yes)

      include "atparam.h"

      include "com_cpl_flags.h"

      include "com_cli_sea.h"
      include "com_var_sea.h"
      
      integer i,j,k
      COMPLEX FIELDS(MX,NX)
      real    sstg(ix,il)


C--   1. Compute climatological fields for initial date

      CALL ATM2SEA (0)

C--   2. Initialize prognostic variables of ocean/ice model
C--      in case of no restart or no coupling

      if (istart.le.0) then
        fields(:,:)=0d0
        fields(3,2)=10
        call grid(fields,sstg,1)
        k=0
        do j=1,nlat
          do i=1,nlon
            k=k+1
            sst_om(k)  = 255d0 + 56*cos(deglat_s(j)/180*3.14)**2 + 20*sin(3*i*2*3.14/real(nlon)) ! SST 
c            sst_om(k) = sstg(i,j)+270
c            write(*,*) sst_om(k)
            tice_om(k) = 273.15     ! sea ice temperature
            sice_om(k) = 0d0     ! sea ice fraction
          enddo
        enddo
c        sst_om(:) = 273d0
c        tice_om(:)= sst_om(:)
c        sice_om(:)=0d0
        
      endif
      
      CALL SEA2ATM (0)

      RETURN
      END


      SUBROUTINE ATM2SEA (jday)
C--
C--   SUBROUTINE ATM2SEA (jday)
C-- 
      include "atparam.h"

      include "com_date.h"
      include "com_cpl_flags.h"

      include "com_cli_sea.h" 
      include "com_var_sea.h"
      include "com_flx_sea.h"

      include "com_cplvar_sea.h"

      real fmasks(ngp)                  ! sea fraction
      equivalence (fmasks,fmask_s)

      real hfyearm(ngp)                 ! annual mean heat flux into the ocean
      equivalence (hfyearm,hfseacl)


C--   1. Interpolate climatological fields and obs. SST anomaly
C--      to actual date

C     Adjust climatological fields over sea ice

c     SST at freezing point
      sstfr = 273.2-1.8

      if (jday.le.0) RETURN

C--   2. Set input variables for mixed-layer/ocean model

      if (icsea.gt.0.or.icice.gt.0) then

        VSEA_INPUT(:,1) = sst_om(:)
        VSEA_INPUT(:,2) = tice_om(:)
        VSEA_INPUT(:,3) = 0d0
        VSEA_INPUT(:,4) = hflux_s(:)
        VSEA_INPUT(:,5) = 0d0
        VSEA_INPUT(:,6) = 283d0
        VSEA_INPUT(:,7) = 273d0
        VSEA_INPUT(:,8) = 0d0

      endif

C--   3. Call message-passing routines to send data (if needed)

      RETURN
      END


      SUBROUTINE SEA2ATM (jday)
C--
C--   SUBROUTINE SEA2ATM (jday)
C-- 
      include "atparam.h"

      include "com_cpl_flags.h"

      include "com_var_sea.h"

      include "com_cplvar_sea.h"

      if (jday.gt.0.and.(icsea.gt.0.or.icice.gt.0)) then

C--   1. Run ocean mixed layer or 
C--      call message-passing routines to receive data from ocean model

         CALL SEA_MODEL 

C--   2. Get updated variables for mixed-layer/ocean model

         sst_om(:)   = VSEA_OUTPUT(:,1)      ! SST
         tice_om(:)  = VSEA_OUTPUT(:,2)      ! sea ice temperature 
         sice_om(:)  = VSEA_OUTPUT(:,3)      ! sea ice fraction

c        sice_om(:)  = sicecl_ob(:)

      endif

C--   3. Compute sea-sfc. anomalies and full fields for atm. model

C     3.1 SST

      sstan_am(:) = 0.

C        Use full ocean model SST
         sst_am(:) = sst_om(:)

C     3.2 Sea ice fraction and temperature

         sice_am(:) = sice_om(:)
         tice_am(:) = tice_om(:)

      ssti_om(:) = sst_om(:)

      RETURN
      END


      SUBROUTINE REST_SEA (imode)
C--
C--   SUBROUTINE REST_SEA (imode)
C--
C--   Purpose : read/write sea variables from/to a restart file
C--   Input :   IMODE = 0 : read model variables from a restart file
C--                   = 1 : write model variables  to a restart file

C-- 
      include "atparam.h"

      include "com_cpl_flags.h"

      include "com_var_sea.h"

      real sst_c(ngp)              ! SST corrected for sea-ice values

      if (imode.eq.0) then

         read (3)  sst_om(:)       ! SST 
         read (3) tice_om(:)       ! sea ice temperature
         read (3) sice_om(:)       ! sea ice fraction

      else

C        write sea/ice model variables from coupled runs,
C        otherwise write fields used by atmospheric model

         sstfr = 273.2-1.8

         if (icsea.gt.0) then
            write (10) sst_om(:) 
         else
            sst_c(:) = max(sst_am(:),sstfr)
            write (10) sst_c(:)
         endif

         if (icice.gt.0) then
            write (10) tice_om(:) 
            write (10) sice_om(:) 
         else
            write (10) tice_am(:)
            write (10) sice_am(:) 
         endif

      endif

      RETURN
      END


      SUBROUTINE OBS_SSTA 
C--
C--   SUBROUTINE OBS_SSTA 
C--
C--   Purpose : update observed SST anomaly array
 
      include "atparam.h"

      include "com_cli_sea.h"

      real*4 r4inp(ix,il)
   
        do j = 1,il
          do i = 1,ix
             sstan3(i,j,1) = sstan3(i,j,2)
             sstan3(i,j,2) = sstan3(i,j,3)
          enddo
        enddo

        read(30,end=100) ((r4inp(i,j),i=1,ix),j=il,1,-1)

        do j = 1,il
          do i = 1,ix
             sstan3(i,j,3) = r4inp(i,j)
          enddo
        enddo 

      CALL FORCHK (bmask_s,sstan3(1,1,3),ix*il,1,-50.,50.,0.)

      RETURN

 100  continue

      print *, ' WARNING: end-of-file reached on SST anomaly file'
      print *, ' SST anomaly will be kept constant'

C--
      RETURN
      END  


