C--
C--   Year, date, time parameters associated with a planet
C--IO s Year, date, time parameters associated with a planet
C--IO s (common TIMEC0)

      MONTHS = 12
      DAYSYR = 365
      DATA DAYSMN/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
     &              0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     &              0,  0,  0,  0,  0,  0,  0,  0 /
      MINSDY = 1440
      SECSDY = 86400
      SECSHR = 3600
      WINTER = 10.

C     Day number of Vernal Equinox, for earth 21 March (30+30+21 in radiat)
      NVE    = 81
C     One year has 360 days for phy_radiat.f = 12 * 30
      DAYRAD = 360
