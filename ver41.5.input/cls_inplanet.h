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
C     SDEP1  = First soil depth for computation of soil water availability
C     IDEP2  = 2nd soil depth factor for computation soil water availability
C     NPOWHD = Power of Laplacian in horizontal diffusion
      SSTFR  = 273.2-1.8
      FRWTR1 = 273.
      FRWTR2 = 273.16
      TGRVWV = 288.
      SDEP1  = 70.
      IDEP2  = 3
      NPOWHD = 4
