C--
C--   Year, date, time parameters associated with a planet
C--IO sx Year, date, time parameters associated with a planet
C--IO sx (common TIMEC0)

      MONTHS = 12
      DAYSYR = 365
      DATA DAYSMN/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
     &              0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &              0,  0,  0,  0,  0,  0,  0,  0 /
      DATA LBLMON /'1jan','1feb','1mar','1apr','1may','1jun',
     &             '1jul','1aug','1sep','1oct','1nov','1dec',
     &             '----','----','----','----','----','----',
     &             '----','----','----','----','----','----',
     &             '----','----','----','----','----','----',
     &             '----','----'/
      MINSDY = 1440
      MINSHR = 60
      SECSDY = 86400
      SECSHR = 3600
      WINTER = 10.

C     For calculation of solar radiation as fn. of latitide
C     Reference: A. Berger, JAS, 35, 2362-2367,1978
C
C     One year has 360 days for phy_radiat.f subroutine solar - 12 * 30
      DAYRAD = 360
C     Day number of Vernal Equinox, for earth 21 March (30+30+21)
      NVE    = 81
C--
C--   Physical parameters associated with a planet
C--IO sx Eccentricity, obliquity, angle between Vernal Equinox and Perihelion
C--IO sx (common PLNTPH)
C     Acceleration due to gravity (m per s^2)
      GRAVIT = 9.81
C     Radius of planet
      REARTH = 6.371E+6
C     pressure vertical velocity [Pa/s]
      OMEGA  = 7.292E-05
C     SEB Physical constant required by the dynamical core
      AKAP   = 2./7.
C     SEB Physical constant required by the dynamical core
      RGAS   = AKAP*1004.
C     Eccentricity (radians)
      ECC    = 0.016724
C     Obliquity (degrees)
      OBLIQ  = 23.446
C     Angle between Vernal Equinox and Perihelion (degrees)
      OMDEG  = 282.04
C     Less accurate obliquity (degrees)
      OBLIQ2 = 23.45
C     Eddy Kinetic energy growth rate maximum
      GRMAX = 0.2/(86400.*2.)
C     Heat capacity of soil per m^2
      HCAPSO = 2.50e+6
C     Heat capacity of ice per m^2
      HCAPIC = 1.93e+6
C     Heat capacity of sea (mixed layer) per m^2
      HCAPSE = 4.18e+6

C     SSTFR  = Sea surface freezing temperature (K)
C              in cpl_sea(_model).f, ppo_dmflux
C     FRWTR1 = Freezing point of 'water' (K) in ini_inbcon.f
C     FRWTR2 = Freezing point of 'water' (K) in phy_convmf.f, phy_lscond.f
C     FRWTR3 = Freezing point of 'water' (K) in phy_shtorh
C     TGRVWV = Temperature used for implicit gravity wave computation
C     TTROP  = tropos:  T = 288 degK at z = 0, constant lapse rate
C     TSTRAT = stratos: T = 216 degK, lapse rate = 0
C     PREF   = Reference pressure for temperature profile (hPa) at z = 0
C     SDEP1  = First soil depth for computation of soil water availability
C     IDEP2  = 2nd soil depth factor for computation soil water availability
C     NPOWHD = Power of Laplacian in horizontal diffusion
C     NPLEVS = Number of standard prseeure levels for post-proc.
C     PLEV   = Standard prseeure levels for post-proc.
C     ESREF  = Reference specific humidity at surface (g/kg)
C     WTRAIR = Ratio water-vapour density to dry air
C     ARMFAC = August-Roche-Magnus factor humidity model
C     ARMC1  = A-R-M constant above freezing humidity model
C     ARMC2  = A-R-M constant below freezing humidity model
C     ARMT1  = A-R-M temperature above freezing humidity model
C     ARMT2  = A-R-M temperature below freezing humidity model
C     BARLPS = Lapse rate for dry air at sea level (K/m)
C     WREF   = Ref. lapse rate to correct temp. exptrapolation - ppo_tminc.f
C     MINST  = Minimum temperature used to calculate mean-sea-level pressure
C     MAXST  = Maximum temperature used to calculate mean-sea-level pressure
      SSTFR  = 273.2-1.8
      FRWTR1 = 273.
      FRWTR2 = 273.16
      TGRVWV = 288.
      TTROP  = 288.
      TSTRAT = 216.
      PREF   = 1013.
      SDEP1  = 70.
      IDEP2  = 3
      NPOWHD = 4
      NPLEVS = 14
      DATA PLEV/ 0.925, 0.850, 0.775, 0.700, 0.600, 0.500, 0.400,
     *           0.300, 0.250, 0.200, 0.150, 0.100, 0.050, 0.030,
     *           0.0,   0.0,   0.0,   0.0,   0.0,   0.0/
      ESREF  = 17.
      WTRAIR = .622
      ARMFAC = 6.108E-3
      ARMC1  = 17.269
      ARMC2  = 21.875
      ARMT1  = 35.86
      ARMT2  =  7.66
      BARLPS = -0.006
      WREF   = 0.7
      MINST  = 255.
      MAXST  = 295.
