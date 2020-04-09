C--:  /TIMEC0/: Time constants for planet (initial. in INTIME)
C--:   WINTER = Difference between day 0 and day of year number of
C--:            winter solstice - +10. days for earth
C--:   MONTHS = Number of months in a year - maximum number 32
C--:   DAYSYR = Number of days in a year
C--:   MINSDY = Number of minutes in a day
C--:   MINSHR = Number of minutes in a day
C--:   SECSDY = Number of seconds in a day
C--:   SECSHR = Number of seconds in an hour
C--:   NVE    = Day of year number of Vernal Equinox
C--:   DAYRAD = Days in a year for phy_radiat.f
C--:   DAYSMN = Number of days in each month (common DAYSMN)
C--:   LBLMON = Labels of 4 characters for each month

      COMMON /TIMEC0/ WINTER, MONTHS, DAYSYR, MINSDY, MINSHR,
     &                SECSDY, SECSHR, NVE, DAYRAD, DAYSMN, LBLMON
C
C--:  /PLNTPH/: Physical parameters associated with a planet
C--:   GRAVIT = Acceleration due to gravity (m per s^2)
C--:   REARTH = Radius of the planet (m)
C--:   OMEGA  = Pressure vertical velocity (Pa/s)
C--:   AKAP   = SEB Physical constant required by the dynamical core
C--:   RGAS   = SEB Physical constant required by the dynamical core
C--:   ECC    = Eccentricity (radians)
C--:   OBLIQ  = Obliquity (degrees)
C--:   OMDEG  = Angle between Vernal Equinox and Perihelion (degrees)
C--:   OBLIQ2 = Less accurate Obliquity (degrees)
C--:   GRMAX  = Eddy Kinetic energy growth rate maximum
C--:   HCAPSO = Heat capacity of soil per m^2
C--:   HCAPIC = Heat capacity of ice per m^2
C--:   HCAPSE = Heat capacity of sea (mixed layer) per m^2
C--:   SSTFR  = Freezing point of 'water' (K) in cpl_sea(_model).f, ppo_dmflux
C--:   FRWTR1 = Freezing point of 'water' (K) in ini_inbcon.f
C--:   FRWTR2 = Freezing point of 'water' (K) in phy_convmf.f, phy_lscond.f
C--:   FRWTR3 = Freezing point of 'water' (K) in phy_shtorh
C--:   TGRVWV = Temperature used for implicit gravity wave computation
C--:   SDEP1  = First soil depth for computation of soil water availability
C--:   IDEP2  = 2nd soil depth factor for computation soil water availability
C--:   NPOWHD = Power of Laplacian in horizontal diffusion
C--:   NPLEVS = Number of standard prseeure levels for post-proc.
C--:   PLEV   = Standard prseeure levels for post-proc.

      COMMON /PLNTPH/ GRAVIT, REARTH, OMEGA, AKAP, RGAS, ECC, OBLIQ,
     $                OMDEG, OBLIQ2, GRMAX, HCAPSO, HCAPIC, HCAPSE,
     $                SSTFR, FRWTR1, FRWTR2, FRWTR3, TGRVWV, PLEV,
     $                SDEP1, IDEP2, NPOWHD, NPLEVS
