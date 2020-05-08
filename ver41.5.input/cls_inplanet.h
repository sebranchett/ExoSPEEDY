C--
C--   Constants for planet parameterisation

C--   Constants for planet calendar (common TIMEC0)

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

      DAYRAD = 360
      NVE    = 81

C--   Physical parameters associated with a planet (common PLNTPH)

      GRAVIT = 9.81
      REARTH = 6.371E+6
      OMEGA  = 7.292E-05
      AKAP   = 2./7.
      RGAS   = AKAP*1004.
      ECC    = 0.016724
      OBLIQ  = 23.446
      OMDEG  = 282.04
      OBLIQ2 = 23.45
      GRMAX  = 0.2/(86400.*2.)

C--   Soil and sea parameters for the planet

      HCAPSO = 2.50e+6
      HCAPIC = 1.93e+6
      HCAPSE = 4.18e+6
      SSTFR  = 273.2-1.8
      FRWTR1 = 273.
      FRWTR2 = 273.16
      FRWTR3 = 273.15
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
      ANOM0  = 20.
