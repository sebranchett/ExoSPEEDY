C--
C--   Length of the integration and time stepping constants (common ISTEPS)
C--IO sx Length of the integration and time stepping constants
C--IO sx Logical flags
C--IO sx 12 months in a year

      NMONTS = 3
      NDAYSL = 0
      NSTEPS = 36

      NSTDIA = 36*5
      NSTPPR = 6
      NSTOUT = -1
      IDOUT  = 0
      NMONRS = 36

      ISEASC = 1
      IYEAR0 = 2000
      IMONT0 = 1

      NSTRAD = 3
      NSTRDF = 0
      INDRDF = 1

      ICLAND = 99
      ICSEA  = 2
      ICICE  = 1
      ISSTAN = 0

      ISSTY0 = 1854
      ISST0  = (IYEAR0-ISSTY0)*MONTHS+IMONT0

C--
C--   Logical flags (common LFLAG1)

      LPPRES = .true.
      LCO2 = .false.
      LTXTO = .true.

