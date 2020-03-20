C--IO sx MAXMON months in a year
C--
C--   /LMASKS/ Land masks
      common /LMASKS/ fmask_l, bmask_l
      real fmask_l(ix,il)                ! fraction of land
      real bmask_l(ix,il)                ! binary land mask
C--
C--   /LCLIM/ Monthly-mean climatological fields over land
      common /LCLIM/ stlmnth, snowdmnth, soilwmnth
      real   stlmnth(ix,il,MAXMON)        ! land surface temperature
      real snowdmnth(ix,il,MAXMON)        ! snow depth (water equiv.)
      real soilwmnth(ix,il,MAXMON)        ! soil water availability
