
      SUBROUTINE INIFFT
C--IO h atparam.h

C     Initialize FFTs

      include "atparam.h"

      COMMON /FFTCOM/ CH(IX),WA(IX),IFAC(IX)

      CALL RFFTI (IX,WA,IFAC)
C      CALL RFFTI (IX,A,IFAC)

      RETURN
      END

*********************************************************************

      SUBROUTINE GRIDX(VARM,VORG,KCOS)
C--IO h atparam.h, com_spectral.h

C     From Fourier coefficients to grid-point data

      include "atparam.h"

      include "com_spectral.h"

      COMMON /FFTCOM/ CH(IX),WA(IX),IFAC(IX)

      REAL VORG(IX,IL), VARM(MX2,IL)
      REAL FVAR(IX)

C     Copy Fourier coefficients into working array

      DO J=1,IL

        FVAR(1)=VARM(1,J)

        DO M=3,MX2
          FVAR(M-1)=VARM(M,J)
        ENDDO
        DO M=MX2,IX
          FVAR(M)=0.0
        ENDDO


C     Inverse FFT

        CALL RFFTB (IX,FVAR,CH,WA,IFAC)
C      CALL RFFTB (IX,FVAR,CH,WA,IFAC)

C     Copy output into grid-point field, scaling by cos(lat) if needed

        IF (KCOS.EQ.1) THEN

          DO I=1,IX
            VORG(I,J)=FVAR(I)
          ENDDO

        ELSE

          DO I=1,IX
            VORG(I,J)=FVAR(I)*COSGR(J)
          ENDDO

        ENDIF

      ENDDO

      RETURN
      END

******************************************************************

      SUBROUTINE SPECX(VORG,VARM)
C--IO h atparam.h

C     From grid-point data to Fourier coefficients

      include "atparam.h"

      COMMON /FFTCOM/ CH(IX),WA(IX),IFAC(IX)

      REAL VORG(IX,IL), VARM(MX2,IL)
      REAL FVAR(IX)

C     Copy grid-point data into working array

      DO J=1,IL

        DO I=1,IX
          FVAR(I)=VORG(I,J)
        ENDDO

C     Direct FFT

        CALL RFFTF (IX,FVAR,CH,WA,IFAC)
C      CALL RFFTF (IX,FVAR,CH,WA,IFAC)

C     Copy output into spectral field, dividing by no. of long.

      SCALE=1./FLOAT(IX)

C     Mean value (a(0))

      VARM(1,J)=FVAR(1)*SCALE
      VARM(2,J)=0.0

      DO M=3,MX2
        VARM(M,J)=FVAR(M-1)*SCALE
      ENDDO

      ENDDO

      RETURN
      END

      include "spe_subfft_fftpack2.f"
