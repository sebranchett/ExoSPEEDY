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
C--:   ECC    = Eccentricity (radians)
C--:   OBLIQ  = Obliquity (degrees)
C--:   OMDEG  = Angle between Vernal Equinox and Perihelion (degrees)
C--:   OBLIQ2 = Less accurate Obliquity (degrees)
C--:   GRMAX  = Eddy Kinetic energy growth rate maximum

      COMMON /PLNTPH/ ECC, OBLIQ, OMDEG, OBLIQ2, GRMAX
