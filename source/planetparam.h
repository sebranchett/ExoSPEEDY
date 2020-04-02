C--IO planet parameters and variable type declarations
      INTEGER     MAXMON
      PARAMETER  (MAXMON=32)
      REAL        WINTER
      INTEGER     MONTHS, DAYSYR, MINSDY, MINSHR,
     &            SECSDY, SECSHR, NVE, DAYRAD
      INTEGER     DAYSMN(MAXMON)
      CHARACTER*4 LBLMON(MAXMON)
      REAL        GRAVIT, REARTH, OMEGA, AKAP, RGAS, ECC, OBLIQ
      REAL        OMDEG, OBLIQ2, GRMAX, HCAPSO, HCAPIC, HCAPSE
