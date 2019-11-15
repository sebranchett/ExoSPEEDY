C--
C--:  /ISTEPS/: Integer time-stepping constants (initial. in INDYNS)
C--:   NMONTS = Integration length in months
C--:   NDAYSL = No. of days in the last month of int. (max=30)
C--:   NSTEPS = No. of time steps in one day
C--:   NSTDIA = Period (no. of steps) for diagnostic print-out
C--:   NSTPPR = Period (no. of steps) for post-processing 
C--:   NSTOUT = Period (no. of steps) for time-mean output
C--:   IDOUT  = daily output flag 
C--:            (0=no, 1=basic (Z500,PREC,MSLP,TEMP0), 2=full)
C--:   NMONRS = Period (no. of months) for restart file update
C--:   ISEASC = Seasonal cycle flag (0=no, 1=yes)
C--    ISTART = Start flag (0: from rest, 1: from restart file)
C--:   IYEAR0 = Year of initial date (4-digit, eg 1900)
C--:   IMONT0 = Month of initial date (1 to 12)

      COMMON /ISTEPS/ NMONTS, NDAYSL, NSTEPS,
     &                NSTDIA, NSTPPR, NSTOUT, IDOUT, NMONRS,
     &                ISEASC, ISTART, IYEAR0, IMONT0

C--
C--:  /ISTFOR/: Integer forcing indices (initial. in INDYNS)
C--:   NSTRAD = Period (no. of steps) for shortwave radiation 
C--:   NSTRDF = Duration of random diabatic forcing ( 0 : no forcing, 
C--:            > 0 : no. of initial steps, < 0 : whole integration)									  
C--:   INDRDF = Initialization index for random diabatic forcing
C--:   ISST0  = record in SST anomaly file corr. to the initial month     

      COMMON /ISTFOR/ NSTRAD, NSTRDF, INDRDF, ISST0

C--
C--   /RSTEPS/: Real time-stepping constants (initial. in INDYNS)
C--    DELT   = Time step in seconds
C--    DELT2  = 2 * time step in seconds
C--    ROB    = Damping factor in Robert time filter
C--    WIL    = Parameter of Williams filter
C--    ALPH   = Coefficient for semi-implicit computations

      COMMON /RSTEPS/ DELT, DELT2, ROB, WIL, ALPH

