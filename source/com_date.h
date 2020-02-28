C--IO sx concpet of 12 months in a year
C--
C--   /DATE1/: date and time variables (updated in NEWDATE)
      common /DATE1/ IYEAR, IMONTH, IDAY, IMONT1, TMONTH, TYEAR, 
     &               RDAY

C--   /DATE2/: calendar set-up (initialized in NEWDATE)
      common /DATE2/ NDAYCAL(MAXMON,2), NDAYTOT
