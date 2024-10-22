C--
C--   Constants for physical parametrization routines:
C--IO sx Constants for physical parametrization routines:
C--IO sx Soil moisture and snow, convection, large-scale condensation
C--IO sx radiation, surface fluxes, vertical diffusion and sh. conv.

C--   Planet physical constants (common PHYCON)

      P0 = 1.e+5
      RD = 287.
      CP = 1004.
      ALHC = 2501.0
cfk      ALHS = 2801.0
      SBC = 5.67e-8

C--   Soil moisture and snow parameters (common SOILMP):

      SWCAP = 0.30
      SWWIL = 0.17
      SD2SC = 60.

C--   Constants for convection (common CNVCON):

      PSMIN  =  0.8
      TRCNV  =  6.
      RHBL   =  0.9
      RHIL   =  0.7
      ENTMAX =  0.5
      SMF    =  0.8

C--   Constants for large-scale condensation (common LSCCON):

      TRLSC  = 8.
      RHLSC  = 0.9
      DRHLSC = 0.1
      RHBLSC = 0.95

C--   Constants for radiation (common RADCON):

      SOLC   = 342.

      ALBSEA = 0.07
      ALBICE = 0.60
      ALBSN  = 0.60

      RHCL1  =  0.30
      RHCL2  =  1.00
      QACL   =  0.20
      WPCL   =  0.2
      PMAXCL = 10.0

      CLSMAX  = 0.60
      CLSMINL = 0.15
      GSE_S0  = 0.25
      GSE_S1  = 0.40

      ALBCL  =  0.43
      ALBCLS =  0.50

      EPSSW  =  0.020
      EPSLW  =  0.05
      EMISFC =  0.98

      ABSDRY =  0.033
      ABSAER =  0.033
      ABSWV1 =  0.022
      ABSWV2 = 15.000

      ABSCL1 =  0.015
      ABSCL2 =  0.15

      ABLWIN =  0.3
      ABLCO2 =  6.0
      RYRCO2  = 1950
      DELCO2 =  0.005
      ABLWV1 =  0.7
      ABLWV2 = 50.0

      ABLCL1 = 12.0
      ABLCL2 =  0.6

C--   Constants for surface fluxes (common SFLCON):

      FWIND0 = 0.95
      FTEMP0 = 1.
      FHUM0  = 0.

      CDL = 2.4e-3
      CDS = 1.0e-3
      CHL = 1.2e-3
      CHS = 0.9e-3

      VGUST  = 5.
      CTDAY  = 1.0e-2
      DTHETA = 3.
      FSTAB  = 0.67
      HDRAG  = 2000.
      FHDRAG = 0.5

      CLAMBDA = 7.
      CLAMBSN = 7. 

C--   Constants for vertical diffusion and sh. conv. (common VDICON):

      TRSHC  =  6.
      TRVDI  = 24.
      TRVDS  =  6.

      REDSHC =  0.5
      RHGRAD =  0.5
      SEGRAD =  0.1
