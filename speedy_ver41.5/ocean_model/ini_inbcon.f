 
      SUBROUTINE INBCON (grav0,radlat)
C--
C--   SUBROUTINE INBCON (grav0,radlat)
C--
C--   Purpose : Read topography and climatological boundary conditions 
C--   Input :   grav0  = gravity accel.
C--             radlat = grid latitudes in radiants
 									
      include "atparam.h"

      include "com_tsteps.h" 
      include "com_cpl_flags.h"
 
      include "com_surfcon.h"    

      include "com_cli_land.h" 
      include "com_cli_sea.h" 

      real   radlat(il)

      real*4 r4inp(ix,il), dummy4
      real*4 veg(ix,il), swl1(ix,il), swl2(ix,il)

      iitest=1

c     set threshold for land-sea mask definition
c     (ie minimum fraction of either land or sea)

      thrsh = 0.1

C--   1. Read topographical fields (orography, land-sea mask)

      if (iitest.ge.1) print*,' read orography' 

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)
      do j = 1,il
        do i = 1,ix
          phi0(i,j) = grav0*r4inp(i,j)*0.00
c *FS*    impose flat surface
c					phi0(i,j) = 0d0
        enddo
      enddo

      call truncg (ntrun,phi0,phis0)
 
      if (iitest.ge.1) print*,' read fractional land-sea mask'  

      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)
      do j = 1,il
        do i = 1,ix
          fmask(i,j) = r4inp(i,j)
        enddo
      enddo

C--   2. Initialize land-sfc boundary conditions

C--   2.1 Fractional and binary land masks

      do j=1,il
       do i=1,ix

         fmask_l(i,j) = fmask(i,j)

         if (fmask_l(i,j).ge.thrsh) then
           bmask_l(i,j) = 1.
           if (fmask(i,j).gt.(1.-thrsh)) fmask_l(i,j) = 1.
         else
           bmask_l(i,j) = 0.
           fmask_l(i,j) = 0.
         endif

c *FS* impose only sea no land

         bmask_l(i,j) = 0.
         fmask_l(i,j) = 0.
         fmask(i,j) = 0d0

         fmask1(i,j) = fmask_l(i,j)

       enddo
      enddo

C--   2.2 Annual-mean surface albedo

      if (iitest.ge.1) print*,' read surface albedo' 
 
      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)
      do j = 1,il
        do i = 1,ix
          alb0(i,j) = 0.01*r4inp(i,j)
          alb0(i,j) = 0.3
        enddo
      enddo


C--   3. Initialize sea-sfc boundary conditions

C--   3.1 Fractional and binary sea masks

      do j=1,il
       do i=1,ix

         fmask_s(i,j) = 1.-fmask(i,j)

         if (fmask_s(i,j).ge.thrsh) then
           bmask_s(i,j) = 1.
           if (fmask_s(i,j).gt.(1.-thrsh)) fmask_s(i,j) = 1.
         else
           bmask_s(i,j) = 0.
           fmask_s(i,j) = 0.
         endif

       enddo
      enddo

C     Grid latitudes for sea-sfc. variables
      rad2deg = 90./asin(1.)
      do j=1,il
         deglat_s(j) = rad2deg*radlat(j)
      enddo

C--   3.2 SST 


      do j = 1,il
        do i = 1,ix
          hfseacl(i,j) = 0.
        enddo
      enddo

C--
      RETURN
      END

      SUBROUTINE FORCHK (FMASK,FIELD,NGP,NF,FMIN,FMAX,FSET)
										
C--   Aux. routine FORCHK: Check consistency of sfc fields with land-sea mask 
C--   and set undefined values to a constant (to avoid over/underflow)

      real fmask(ngp), field(ngp,nf)

      do jf = 1,nf

        nfault=0

        do jgp = 1,ngp
          if (fmask(jgp).gt.0.0) then
            if (field(jgp,jf).lt.fmin.or.field(jgp,jf).gt.fmax)
     *          nfault = nfault+1
          else
            field(jgp,jf) = fset
          endif
        enddo

        print *, ' field: ', jf, '   no. of faulty points:', nfault

      enddo

      print *, ' undefined values set to', fset

      RETURN
      END


      SUBROUTINE TRUNCG (ITR,FG1,FG2)

C--   SUBROUTINE TRUNCG (ITR,FG1,FG2)
C--   Purpose : compute a spectrally-filtered grid-point field
C--   Input   : ITR : spectral truncation (triangular)
C--           : FG1 : original grid-point field
C--   Output  : FG2 : filtered grid-point field

      include "atparam.h"

      REAL FG1 (IX,IL), FG2(IX,IL)
      COMPLEX FSP(MX,NX), ZERO 

      ZERO = (0.,0.)

      CALL SPEC (FG1,FSP)

      DO N=1,NX
        DO M=1,MX
          ITWN=ISC*(M-1)+N-1
          IF (ITWN.GT.ITR) FSP(M,N)=ZERO
        ENDDO
      ENDDO

      CALL GRID (FSP,FG2,1)

      RETURN
      END
