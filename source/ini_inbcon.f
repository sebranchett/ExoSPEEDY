 
      SUBROUTINE INBCON (grav0,radlat)
C--
C--   SUBROUTINE INBCON (grav0,radlat)
C--
C--   Purpose : Read topography and climatological boundary conditions 
C--   Input :   grav0  = gravity accel.
C--             radlat = grid latitudes in radiants
C--IO h atparam.h, com_tsteps.h, com_cpl_flags.h
C--IO h com_surfcon.h, com_cli_land.h, com_cli_sea.h
C--IO h com_planet.h
C--IO sx 12 months in a year
C--IO sx sdep1 = 70. and idep2 = 3, soil depths
C--IO sx 273., temperature of freezing water? FRWTR
C--IO r read topographical fields from unit 20 - orography
C--IO r read topographical fields from unit 20 - land-sea mask
C--IO r read topographical fields from unit 20 - albedo
C--IO r read land-surface temp. from unit 23
C--IO r read snow depth from unit 24
C--IO r read vegetation fraction from unit 25
C--IO r read soil moisture from unit 26
C--IO r read 'SST' (sea surface temperature?) from unit 21
C--IO r read sea ice fraction from unit 22
C--IO r read SST anomalies initial and prec./following months from unit 30
C--IO r read Annual-mean heat flux into sea-surface from unit 31
C--IO r read ocean model SST bias from unit 32 - commented out
C--IO w write correction for model-to-actual topography to unit 18

      include "atparam.h"
      include "planetparam.h"

      include "com_planet.h"
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
          phi0(i,j) = grav0*r4inp(i,j)
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

         fmask1(i,j) = fmask_l(i,j)

       enddo
      enddo

C--   2.2 Annual-mean surface albedo

      if (iitest.ge.1) print*,' read surface albedo' 
 
      read (20) ((r4inp(i,j),i=1,ix),j=il,1,-1)
      do j = 1,il
        do i = 1,ix
          alb0(i,j) = 0.01*r4inp(i,j)
        enddo
      enddo

C--   2.3 Land-surface temp.

      if (iitest.ge.1) print*,' reading land-surface temp.'
  
      do it = 1,months
        read (23) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            stlmnth(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.eq.1) print*,' checking land-surface temp.'

      CALL FORCHK (bmask_l,stlmnth,ix*il,months,0.,400.,FRWTR)

C     Correction for model-to-actual topography
      do it = 1,months

        call ftland (stlmnth(1,1,it),phi0,phis0,bmask_l)

        if (iitest.gt.1) then
          do j = 1,il
            do i = 1,ix
              r4inp(i,j) = stlmnth(i,j,it)
            enddo
          enddo
          write (18) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        endif

      enddo

C--   2.4 Snow depth

      if (iitest.ge.1) print*,' reading snow depth'  

      do it = 1,months
        read (24) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            snowdmnth(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking snow depth'

      CALL FORCHK (bmask_l,snowdmnth,ix*il,months,0.,20000.,0.)

C--  2.5 Read soil moisture and compute soil water availability 
C--      using vegetation fraction

      if (iitest.ge.1) print*,' reading soil moisture'  

c     read vegetation fraction (in %)
      read (25) ((veg(i,j),i=1,ix),j=il,1,-1)

      do j = 1,il
        do i = 1,ix
          veg(i,j)=max(0.,0.01*veg(i,j))
        enddo
      enddo

      sdep2 = idep2*sdep1

      swwil2= sdep2*swwil
      rsw   = 1./(sdep1*swcap+sdep2*(swcap-swwil))

      do it = 1,months

        read (26) ((swl1(i,j),i=1,ix),j=il,1,-1)
        read (26) ((swl2(i,j),i=1,ix),j=il,1,-1)
        read (26) dummy4

        do j = 1,il
          do i = 1,ix
            swroot = idep2*swl2(i,j)
            soilwmnth(i,j,it) = min(1.,rsw*(swl1(i,j)+
     &                        veg(i,j)*max(0.,swroot-swwil2)))
          enddo
        enddo

      enddo

      if (iitest.ge.1) print*,' checking soil moisture'

      CALL FORCHK (bmask_l,soilwmnth,ix*il,months,0.,10.,0.)


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

      if (iitest.ge.1) print*,' reading sst' 

      do it = 1,months
        read (21) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            sstmnth(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking sst'

      CALL FORCHK (bmask_s,sstmnth,ix*il,months,100.,400.,FRWTR)

c     3.2 Sea ice fraction

      if (iitest.ge.1) print*,' reading sea ice'  

      do it = 1,months
        read (22) ((r4inp(i,j),i=1,ix),j=il,1,-1)
        do j = 1,il
          do i = 1,ix
            sicemnth(i,j,it) = r4inp(i,j)
          enddo
        enddo
      enddo

      if (iitest.ge.1) print*,' checking sea ice'

      CALL FORCHK (bmask_s,sicemnth,ix*il,months,0.,1.,0.)

C--   3.3 SST anomalies for initial and prec./following months

      if (isstan.gt.0) then

        if (iitest.ge.1) print*,' reading sst anomalies' 

        do jrec=1,isst0-2
          read (30) dummy4
        enddo

        do it=1,3

c         NB If isst0 = 1, SST_an(it=1) is set = SST_an(it=2)  
          if (it.ne.2.or.isst0.gt.1)
     &       read (30) ((r4inp(i,j),i=1,ix),j=il,1,-1)

          do j = 1,il
            do i = 1,ix
              sstan3(i,j,it) = r4inp(i,j)
            enddo
          enddo

        enddo

        if (iitest.ge.1) print*,' checking sst anomalies'

        CALL FORCHK (bmask_s,sstan3,ix*il,3,-50.,50.,0.)

      endif

C--   3.4. Annual-mean heat flux into sea-surface

      do j = 1,il
        do i = 1,ix
          hfseacl(i,j) = 0.
        enddo
      enddo

      if (icsea.ge.1) then

        if (iitest.ge.1) print*,' reading sfc heat fluxes' 

        irecl = 4*ix*il
        irec = 0

        open ( unit=31, file='fort.31', status='old', 
     &         form='unformatted', access='direct', recl=irecl )

        do it = 1,months

          irec=irec+2
          read (31,rec=irec) r4inp

          do j = 1,il
            do i = 1,ix
              hfseacl(i,j) = hfseacl(i,j)+r4inp(i,j)
            enddo
          enddo  

        enddo

        do j = 1,il
          do i = 1,ix
            if (bmask_s(i,j).gt.0.) then
                hfseacl(i,j) = hfseacl(i,j)/(1.*months*fmask_s(i,j))
            else
                hfseacl(i,j) = 0.
            endif
          enddo
        enddo  

        if (iitest.ge.1) print*,' checking sfc heat fluxes'

        CALL FORCHK (bmask_s,hfseacl,ix*il,1,-1000.,1000.,0.)

      endif

C--   3.5. Ocean model SST climatology:
C--        defined by adding SST model bias to obs. climatology
C--        (bias may be defined in a different period from climatology)

      if (icsea.ge.3) then

        if (iitest.ge.1) print*,' reading ocean model SST bias' 

c        irecl = 4*ix*il
c        irec = 0

c        open ( unit=32, file='fort.32', status='old', 
c     &         form='unformatted', access='direct', recl=irecl )

        do it = 1,months

c          irec=irec+1
c          read (32,rec=irec) r4inp
cfk           read (32) r4inp

          do j = 1,il
            do i = 1,ix
cfk              sstommnth(i,j,it) = sstmnth(i,j,it)+r4inp(i,j)
              sstommnth(i,j,it) = sstmnth(i,j,it)
            enddo
          enddo  

        enddo

        if (iitest.ge.1) print*,' checking ocean model SST'

        CALL FORCHK (bmask_s,sstommnth,ix*il,months,100.,400.,FRWTR)

      endif
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

      SUBROUTINE FTLAND (STL,PHI0,PHIS0,FMASKL)
C--IO h atparam.h, atparam1.h, com_dyncon0.h, com_dyncon1.h
C--IO h planetparam.h, com_planet.h

      include "atparam.h"
      include "atparam1.h"
      include "planetparam.h"

      include "com_planet.h"
      include "com_dyncon0.h" 
      include "com_dyncon1.h" 

      REAL STL(NLON,NLAT), PHI0(NLON,NLAT), PHIS0(NLON,NLAT),
     &     FMASKL(NLON,NLAT)

      REAL STL2(NLON,NLAT)

      NL8 = NLAT/8
      GAM = 0.001*GAMMA/GRAVIT

      NLAT1 = 1
      NLAT2 = NL8

      DO JBAND=1,8

        SUMT=0.
        SUMW=0.

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           STL(I,J)=STL(I,J)+GAM*PHI0(I,J)
           SUMT=SUMT+GCOS(J)*FMASKL(I,J)*STL(I,J)
           SUMW=SUMW+GCOS(J)*FMASKL(I,J)
         ENDDO
        ENDDO

        SUMT=SUMT/SUMW

        DO J=NLAT1,NLAT2
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=SUMT
         ENDDO
        ENDDO
  
        NLAT1=NLAT1+NL8
        NLAT2=NLAT2+NL8

      ENDDO

      ITR=7
      IDTR=(NTRUN-6)/3

      DO JFIL=1,4

        CALL TRUNCG (ITR,STL,STL2)

        DO J=1,NLAT
         DO I=1,NLON
           IF (FMASKL(I,J).EQ.0.) STL(I,J)=STL2(I,J)
         ENDDO
        ENDDO

        ITR=MIN(ITR+IDTR,NTRUN)

      ENDDO

      CALL TRUNCG (ITR,STL,STL2)

      DO J=1,NLAT
       DO I=1,NLON
         STL(I,J)=STL2(I,J)-GAM*PHIS0(I,J)
       ENDDO
      ENDDO       

      RETURN
      END

      SUBROUTINE TRUNCG (ITR,FG1,FG2)

C--   SUBROUTINE TRUNCG (ITR,FG1,FG2)
C--   Purpose : compute a spectrally-filtered grid-point field
C--   Input   : ITR : spectral truncation (triangular)
C--           : FG1 : original grid-point field
C--   Output  : FG2 : filtered grid-point field

C--IO h atparam.h
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
