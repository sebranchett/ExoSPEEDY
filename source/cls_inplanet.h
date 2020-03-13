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
C--IO s Eccentricity, obliquity, angle between Vernal Equinox and Perihelion
C--IO s (common PLNTPH)
C     Eccentricity (radians)
      ECC    = 0.016724
C     Obliquity (degrees)
      OBLIQ  = 23.446
C     Angle between Vernal Equinox and Perihelion (degrees)
      OMDEG  = 282.04
C     Less accurate obliquity (degrees)
      OBLIQ2 = 23.45
