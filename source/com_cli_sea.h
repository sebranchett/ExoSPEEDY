
C--IO sx 12 months in a year
C--
C--   /SMASKS/ Sea masks
      common /SMASKS/ fmask_s, bmask_s, deglat_s
      real fmask_s(ix,il)                ! fraction of sea
      real bmask_s(ix,il)                ! binary sea mask
      real deglat_s(il)                  ! grid latitudes 
C--
C--   /SCLIM/ Monthly-mean climatological fields over sea
      common /SCLIM/ sstmnth, sicemnth
      real  sstmnth(ix,il,MAXMON)        ! sea/ice surface temperature
      real sicemnth(ix,il,MAXMON)        ! sea ice fraction
C--
C--   /SSTANOM/ SST anomaly fields
      common /SSTANOM/ sstan3
      real sstan3(ix,il,3)              ! sst anomaly in 3 consecutive months
C--
C--   /SMCLIM/ Climatological fields from model output
      common /SMCLIM/ hfseacl, sstommnth
      real hfseacl(ix,il)               ! annual-mean heat flux into sea sfc.
      real sstommnth(ix,il,MAXMON)      ! ocean model SST climatology
