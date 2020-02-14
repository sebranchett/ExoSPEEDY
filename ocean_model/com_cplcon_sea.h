
      real  rhcaps(nlon,nlat)           ! 1./heat_capacity (sea)
      real  rhcapi(nlon,nlat)           ! 1./heat_capacity (ice)
      real   cdsea(nlon,nlat)           ! 1./dissip_time (sea)
      real   cdice(nlon,nlat)           ! 1./dissip_time (ice)
      real   beta                       ! Heat flux coef. at sea/ice int.
      real  dmpsea(mx,nx)               ! spectral diffusion coefficients
      real  dts                         ! timestep in seconds
      real  rlatdeg(nlat)               ! latitudes in degrees

C--
C--   /SEA_MC/ Constant parameters and fields in sea/ice model
      common /SEA_MC/rhcaps,rhcapi,cdsea,cdice,beta,dmpsea,dts,rlatdeg

