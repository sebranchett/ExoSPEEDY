C--
C--   /DATE1/: date and time variables (updated in NEWDATE)
      common /DATE1/ IYEAR, IMONTH, IDAY, IMONT1, TMONTH, TYEAR, 
     &               RDAY

C--   /DATE2/: calendar set-up (initialized in NEWDATE)
      common /DATE2/ NDAYCAL(12,2), NDAYTOT
