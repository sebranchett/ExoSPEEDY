

      SUBROUTINE SEA_MODEL_INIT (fmask_s,rlat) 
C--
C--   SUBROUTINE SEA_MODEL_INIT (fmask_s,rlat)
C--
C--   Purpose : Initialization of sea model
C--   Initialized common blocks: SEA_MC
C--	
C--IO h atparam.h, atparam1.h, com_cplcon_sea.h, cls_insea.h, com_hdifcon.h
C--IO h cls_indyns.h, planetparam.h, com_planet.h
C--IO sx depth_ml for high-latitutde depth, dept0_ml for minimum depth
C--IO sx depth_ice for high-latitude depth, dept0_ice for minimum depth
C--IO sx tdsst for sea-surface dissipation, tdice for sea-ice dissipation
C--IO sx beta for heat flux coefficient
C--IO sx l_globe for global domain logical
C--IO sx l_northe, l_natlan, l_npacif for regional domain logicals
C--IO sx l_tropic, l_indian for regional domain logicals
C--IO sx heat capacity of mixed-l. = 4.18e+6 used for hcaps
C--IO sx heat capacity of sea-ice = 1.93e+6 used for hcapi
C--IO sx 86400 for reciprocal heat capacities = seconds in a day?

      include "atparam.h"
      include "atparam1.h"
      include "planetparam.h"
      include "com_planet.h"
      include "com_hdifcon.h"
      include "com_cplcon_sea.h"

C     Input variables

      real fmask_s(nlon,nlat)            ! sea mask (fraction of sea)
      real rlat(nlat)                    ! latitudes in degrees

C     Auxiliary variables

      real    dmask(nlon,nlat)              ! domain mask
      logical l_globe,
     &        l_northe, l_natlan, l_npacif, 
     &        l_tropic, l_indian            ! domain flags

      real hcaps(nlat)                      ! heat capacity of mixed-l.
      real hcapi(nlat)                      ! heat capacity of sea-ice

      include "cls_indyns.h"

C--  
C--   1. Set geographical domain, heat capacities and dissipation times
C--      for sea (mixed layer) and sea-ice 
C     Model parameters (default values)
C     time step in seconds

      dts = REAL(SECSDY)

      include "cls_insea.h"

C     Heat capacities per m^2 (depth*heat_cap/m^3)

      crad=asin(1.)/90.
      do j=1,nlat
        coslat   = cos(crad*rlat(j))
        rlatdeg(j) = rlat(j)
        hcaps(j) = HCAPSE*(depth_ml +(dept0_ml -depth_ml) *coslat**3)
        hcapi(j) = HCAPIC*(depth_ice+(dept0_ice-depth_ice)*coslat**2)
      enddo
C--
C--   3. Compute constant parameters and fields

C     Set domain mask

      if (l_globe) then
        dmask(:,:) = 1.
      else
        dmask(:,:) = 0.
        if (l_northe) call SEA_DOMAIN ('northe',rlat,dmask)
cfk        if (l_arctic) call SEA_DOMAIN ('arctic',rlat,dmask)
        if (l_natlan) call SEA_DOMAIN ('natlan',rlat,dmask)
        if (l_npacif) call SEA_DOMAIN ('npacif',rlat,dmask)
        if (l_tropic) call SEA_DOMAIN ('tropic',rlat,dmask)
        if (l_indian) call SEA_DOMAIN ('indian',rlat,dmask)
      endif

C     Set heat capacity and dissipation time over selected domain

      do j=1,nlat
         rhcaps(:,j) = 1./hcaps(j)
         rhcapi(:,j) = 1./hcapi(j)
      enddo

      cdsea(:,:) = dmask(:,:)*tdsst/(1.+dmask(:,:)*tdsst)
      cdice(:,:) = dmask(:,:)*tdice/(1.+dmask(:,:)*tdice)

      dmpsea(:,:)=dmps(:,:)*THDS/(24d0*tdsst)

      return
      end

      SUBROUTINE SEA_MODEL 
C--
C--   SUBROUTINE SEA_MODEL
C--
C--   Purpose : Integrate slab ocean and sea-ice models for one day
C--IO h atparam.h, com_cplcon_sea.h, com_cplvar_sea.h
C--IO h planetparam.h, com_planet.h
C--IO sx SST at freezing point = 273.2-1.8
C--IO sx anom0 = 20. for non-linear damping coefficient

      include "atparam.h"
      include "planetparam.h"

c      real vsea_input(nlon,nlat,8), vsea_output(nlon,nlat,3)

      include "com_planet.h"
      include "com_cplcon_sea.h"

      include "com_cplvar_sea.h"

C     Input variables:

      real  sst0(nlon,nlat)     ! SST at initial time
      real tice0(nlon,nlat)     ! sea ice temp. at initial time
      real sice0(nlon,nlat)     ! sea ice fraction at initial time
      real hfsea(nlon,nlat)     ! sea+ice  sfc. heat flux between t0 and t1
      real hfice(nlon,nlat)     ! ice-only sfc. heat flux between t0 and t1

      real  sstcl1(nlon,nlat)   ! clim. SST at final time 
      real ticecl1(nlon,nlat)   ! clim. sea ice temp. at final time
      real hfseacl(nlon,nlat)   ! clim. heat flux due to advection/upwelling

      equivalence    (sst0,vsea_input(1,1))
      equivalence   (tice0,vsea_input(1,2))
      equivalence   (sice0,vsea_input(1,3))
      equivalence   (hfsea,vsea_input(1,4))
      equivalence   (hfice,vsea_input(1,5))
      equivalence  (sstcl1,vsea_input(1,6))
      equivalence (ticecl1,vsea_input(1,7))
      equivalence (hfseacl,vsea_input(1,8))

C     Output variables
 
      real  sst1(nlon,nlat)     ! SST at final time
      real  sstd(nlon,nlat)     ! diffusive SST change
      real tice1(nlon,nlat)     ! sea ice temp. at final time 
      real sice1(nlon,nlat)     ! sea ice fraction at final time 

      equivalence  (sst1,vsea_output(1,1))
      equivalence (tice1,vsea_output(1,2))
      equivalence (sice1,vsea_output(1,3))

C     Auxiliary variables

      real hflux(nlon,nlat)   ! net sfc. heat flux
      real tanom(nlon,nlat)   ! sfc. temperature anomaly
      real  cdis(nlon,nlat)   ! dissipation ceofficient

c      beta = 1.               ! heat flux coef. at sea-ice bottom
C--
C--   1. Ocean mixed layer

C     Net heat flux
      hflux(:,:) = hfsea(:,:)

      sstd(:,:) = rhcaps(:,:)*hflux(:,:)

      call HORDIFSEA(sst0,sstd)

C     Time evoloution of temp. anomaly

      sst1(:,:) = sst0 + sstd(:,:)*dts

      crad=asin(1.)/90.
      sum=0d0
      asum=0d0
      do j=1,nlat
        coslat   = cos(crad*rlatdeg(j))
        sum=sum+sst1(1,j)*coslat
        asum=asum+coslat
      enddo

      write(*,*) 'mean ocean T ',sum/asum - FRWTR,asum,sum

C--
C--   2. Sea-ice slab model


C     Full ice temperature at final time
      tice1(:,:) = FRWTR

C     Persistence of sea ice fraction
      sice1(:,:) = 0d0

      return
      end

      SUBROUTINE HORDIFSEA (FIELD,FIELDT)
C--
C--   Aux. subr. HORDIF (NLEV,FIELD,FDT,DMP,DMP1)
C--   Purpose : Add horizontal diffusion tendency of FIELD
C--             to grid tendency FIELDT
C--             using damping coefficients DMPSEA
C--
      include "atparam.h"
      include "com_cplcon_sea.h"

      COMPLEX FIELDS(MX,NX), FIELDTS(MX,NX)

      REAL FIELD (IX,IL), FIELDT (IX,IL)

      INTEGER M,N


      CALL SPEC (FIELDT,FIELDTS)
      CALL SPEC (FIELD,FIELDS)

      FIELDTS=(FIELDTS-DMPSEA*FIELDS)/(1+DMPSEA*DTS)

      CALL GRID   (FIELDTS, FIELDT,  1)


      RETURN
      END

      SUBROUTINE SEA_DOMAIN (cdomain,rlat,dmask) 
C--
C--   SUBROUTINE SEA_DOMAIN (cdomain,rlat,dmask)
C--
C--   Purpose : Definition of ocean domains
C--	
C--IO h atparam.h
C--IO sx lat/lon masks set dependent on earth's sea domains

      include "atparam.h"

C     Input variables

      character*6 cdomain           ! domain name
      real rlat(nlat)               ! latitudes in degrees

C     Output variables (initialized by calling routine) 

      real dmask(nlon,nlat)         ! domain mask 

      print *, 'sea domain : ', cdomain

      dlon = 360./float(nlon)
                                   
      if (cdomain.eq.'northe') then
        do j=1,nlat
          if (rlat(j).gt.20.0) dmask(:,j) = 1.
        enddo
      endif

      if (cdomain.eq.'natlan') then
         do j=1,nlat
           if (rlat(j).gt.20.0.and.rlat(j).lt.80.0) then
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.lt.45.0.or.rlon.gt.260.0) dmask(i,j) = 1.
             enddo
           endif
         enddo
      endif

      if (cdomain.eq.'npacif') then
         do j=1,nlat
           if (rlat(j).gt.20.0.and.rlat(j).lt.65.0) then
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.gt.120.0.and.rlon.lt.260.0) dmask(i,j) = 1.
             enddo
           endif
         enddo
      endif

      if (cdomain.eq.'tropic') then
        do j=1,nlat
          if (rlat(j).gt.-30.0.and.rlat(j).lt.30.0) dmask(:,j) = 1.
        enddo
      endif

      if (cdomain.eq.'indian') then
         do j=1,nlat
           if (rlat(j).gt.-30.0.and.rlat(j).lt.30.0) then
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.gt.30.0.and.rlon.lt.120.0) dmask(i,j) = 1.
             enddo
           endif
         enddo
      endif

      if (cdomain.eq.'elnino') then
         do j=1,nlat
           arlat = abs(rlat(j))
           if (arlat.lt.25.0) then
             wlat = 1.
             if (arlat.gt.15.0) wlat = (0.1*(25.-arlat))**2
             rlonw = 300.-2*max(rlat(j),0.)
             do i=1,nlon
               rlon = (i-1)*dlon
               if (rlon.gt.165.0.and.rlon.lt.rlonw) then
                  dmask(i,j) = wlat
               else if (rlon.gt.155.0.and.rlon.lt.165.0) then 
                  dmask(i,j) = wlat*0.1*(rlon-155.)
               endif
             enddo
           endif
         enddo
      endif

c      do j=1,nlat
c       print *, 'lat = ',rlat(j),' sea model domain  = ',dmask(nlon/2,j)
c      enddo

      return
      end

