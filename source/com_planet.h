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
C--:   AKAP   = Poisson constant for dry air at constant pressure
C--             (ratio of gas constant to specific heat)
C--:   RGAS   = Dry gas constant of the atmosphere (J K^-1 kg^-1)
C--:   ECC    = Eccentricity (radians)
C--:   OBLIQ  = Obliquity (degrees)
C--:   OMDEG  = Angle between Vernal Equinox and Perihelion (degrees)
C--:   GRMAX  = Eddy Kinetic energy growth rate maximum
C--:   HCAPSO = Heat capacity of soil per m^2
C--:   HCAPIC = Heat capacity of ice per m^2
C--:   HCAPSE = Heat capacity of sea (mixed layer) per m^2
C--:   SSTFR  = Freezing point of 'water' (K) in cpl_sea(_model).f, ppo_dmflux
C--:   FRWTR  = Freezing point of 'water' (K)
C--:   TGRVWV = Temperature used for implicit gravity wave computation
C--:   TTROP  = tropos:  T = 288 degK at z = 0, constant lapse rate
C--:   TSTRAT = stratos: T = 216 degK, lapse rate = 0
C--:   PREF   = Reference pressure for temperature profile (hPa) at z = 0
C--:   ESREF  = Reference specific humidity at surface (g/kg)
C--:   WTRAIR = Density ratio water-vapour to dry air
C--:   ARMFAC = 6.108E-3 - August-Roche-Magnus factor humidity model
C--:   ARMC1  = 17.269 - A-R-M constant above freezing humidity model
C--:   ARMC2  = 21.875 - A-R-M constant below freezing humidity model
C--:   ARMT1  = 35.86 - A-R-M temperature above freezing humidity model
C--:   ARMT2  =  7.66 - A-R-M temperature below freezing humidity model
C--:   SDEP1  = First soil depth for computation of soil water availability
C--:   IDEP2  = 2nd soil depth factor for computation soil water availability
C--:   NPOWHD = Power of Laplacian in horizontal diffusion
C--:   NPLEVS = Number of standard prseeure levels for post-proc.
C--:   PLEV   = Standard prseeure levels for post-proc.
C--:   BARLPS = Lapse rate for dry air at sea level (K/m)
C--:   WREF   = Ref. lapse rate to correct temp. exptrapolation - ppo_tminc.f
C--:   MINST  = Minimum temperature used to calculate mean-sea-level pressure
C--:   MAXST  = Maximum temperature used to calculate mean-sea-level pressure
C--:   ANOM0  = Non-linear damping coefficient for anomolies

      COMMON /PLNTPH/ GRAVIT, REARTH, OMEGA, AKAP, RGAS, ECC, OBLIQ,
     $                OMDEG, GRMAX, HCAPSO, HCAPIC, HCAPSE,
     $                SSTFR, FRWTR, TGRVWV, PLEV,
     $                TTROP, TSTRAT, PREF, ESREF, WTRAIR,
     $                ARMFAC, ARMC1, ARMC2, ARMT1, ARMT2, BARLPS, WREF,
     $                MINST, MAXST, ANOM0, SDEP1, IDEP2, NPOWHD, NPLEVS
