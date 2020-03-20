C--IO planet parameters and variable type declarations
      INTEGER     MAXMON
      PARAMETER  (MAXMON=32)
      REAL        WINTER
      INTEGER     MONTHS, DAYSYR, MINSDY,
     &            SECSDY, SECSHR, NVE, DAYRAD
      INTEGER     DAYSMN(MAXMON)
      CHARACTER*4 LBLMON(MAXMON)
      REAL        ECC, OBLIQ, OMDEG, OBLIQ2, GRMAX
