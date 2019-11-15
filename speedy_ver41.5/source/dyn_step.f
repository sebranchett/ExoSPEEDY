      SUBROUTINE STEP (J1,J2,DT,ALPH,ROB,WIL)
C--
C--   SUBROUTINE STEP (J1,J2,DT,ALPH,ROB,WIL)
C--
C--   Purpose: perform one time step starting from F(1) and F(2) 
C--            and using the following scheme:
C--
C--   Fnew = F(1) + DT * [ T_dyn(F(J2)) + T_phy(F(1)) ]
C--   F(1) = (1-2*eps)*F(J1) + eps*[F(1)+Fnew]
C--   F(2) = Fnew
C--
C--   Input: 
C--   If J1=1, J2=1 : forward time step (eps=0)
C--   If J1=1, J2=2 : initial leapfrog time step (eps=0)
C--   If J1=2, J2=2 : leapfrog time step with time filter (eps=ROB)
C--   DT = time step (if DT < or = 0, tendencies are computed but 
C--                   no time stepping is performed)
C--   ALPH = 0   : forward step for gravity wave terms
C--   ALPH = 1   : backward implicit step for g.w.
C--   ALPH = 0.5 : centered implicit step for g.w.
C--   ROB  = Robert filter coefficient
C--   WIL  = Williams filter coefficient
C-- 
C--   Modified common blocks : DYNSP1, DYNSP2
C--
      include "atparam.h"
      include "atparam1.h"

      include "com_dyncon0.h"
      include "com_hdifcon.h"

      include "com_dynvar.h"

      COMPLEX VORDT(MX,NX,KX), DIVDT(MX,NX,KX), TDT(MX,NX,KX),
     *        PSDT(MX,NX), TRDT(MX,NX,KX,NTR)

      COMPLEX CTMP(MX,NX,KX)

      iitest=0
      if(iitest.eq.1) print*, ' inside step'

C--   1. Computation of grid-point tendencies
C         (converted to spectral at the end of GRTEND)

      if (iitest.eq.1) print*,' call grtend'
      CALL GRTEND (VORDT,DIVDT,TDT,PSDT,TRDT,1,J2)

C--   2. Computation of spectral tendencies

      IF (ALPH.EQ.0.) THEN

        if (iitest.eq.1) print*,' call sptend'
        CALL SPTEND (DIVDT,TDT,PSDT,J2)

      ELSE

        if (iitest.eq.1) print*,' call sptend'
        CALL SPTEND (DIVDT,TDT,PSDT,1)

c       implicit correction 
        if (iitest.eq.1) print*,' call implic'
        CALL IMPLIC (DIVDT,TDT,PSDT)

      ENDIF

C--   3. Horizontal diffusion

      if (iitest.eq.1) print*, ' biharmonic damping '

C     3.1 Diffusion of wind and temperature
 
      CALL HORDIF (KX,VOR,VORDT,DMP, DMP1)
      CALL HORDIF (KX,DIV,DIVDT,DMPD,DMP1D)

      DO K=1,KX
        DO M=1,MX
          DO N=1,NX 
            CTMP(M,N,K)=T(M,N,K,1)
     &                 +TCORH(M,N)*TCORV(K)
          ENDDO
        ENDDO
      ENDDO

      CALL HORDIF (KX,CTMP,TDT,DMP,DMP1)

C     3.2 Stratospheric diffusion and zonal wind damping

      SDRAG=1./(TDRS*3600.)
      DO N=1,NX
        VORDT(1,N,1)=VORDT(1,N,1)-SDRAG*VOR(1,N,1,1)
        DIVDT(1,N,1)=DIVDT(1,N,1)-SDRAG*DIV(1,N,1,1)
      ENDDO

      CALL HORDIF (1,VOR, VORDT,DMPS,DMP1S)
      CALL HORDIF (1,DIV, DIVDT,DMPS,DMP1S)
      CALL HORDIF (1,CTMP,TDT,  DMPS,DMP1S)

C     3.3 Chech for eddy kinetic energy growth rate 
    
C      CALL CGRATE (VOR,DIV,VORDT,DIVDT)

C     3.4 Diffusion of tracers

      DO K=1,KX
        DO M=1,MX
          DO N=1,NX
           CTMP(M,N,K)=TR(M,N,K,1,1)
     &                 +QCORH(M,N)*QCORV(K)
          ENDDO
        ENDDO
      ENDDO

      CALL HORDIF (KX,CTMP,TRDT,DMPD,DMP1D)

      IF (NTR.GT.1) THEN
        DO ITR=2,NTR
          CALL HORDIF (KX,TR(1,1,1,1,ITR),TRDT(1,1,1,ITR),
     &                 DMP,DMP1)
        ENDDO
      ENDIF

C--   4. Time integration with Robert filter

      IF (DT.LE.0.) RETURN

      if (iitest.eq.1) print*,' time integration'

      IF (J1.EQ.1) THEN
        EPS=0.
      ELSE
        EPS=ROB
      ENDIF

      CALL TIMINT (J1,DT,EPS,WIL,1,PS,PSDT)

      CALL TIMINT (J1,DT,EPS,WIL,KX,VOR,VORDT)
      CALL TIMINT (J1,DT,EPS,WIL,KX,DIV,DIVDT)
      CALL TIMINT (J1,DT,EPS,WIL,KX,T,  TDT)

      DO ITR=1,NTR
        CALL TIMINT (J1,DT,EPS,WIL,KX,
     &               TR(1,1,1,1,ITR),TRDT(1,1,1,ITR))
      ENDDO

C--
      RETURN
      END   

      SUBROUTINE HORDIF (NLEV,FIELD,FDT,DMP,DMP1)
C--
C--   Aux. subr. HORDIF (NLEV,FIELD,FDT,DMP,DMP1)
C--   Purpose : Add horizontal diffusion tendency of FIELD 
C--             to spectral tendency FDT at NLEV levels
C--             using damping coefficients DMP and DMP1
C--
      include "atparam.h"

      COMPLEX FIELD(MXNX,NLEV), FDT(MXNX,NLEV)
      REAL    DMP(MXNX), DMP1(MXNX)

      DO K=1,NLEV
        DO M=1,MXNX
          FDT(M,K)=(FDT(M,K)-DMP(M)*FIELD(M,K))*DMP1(M)
        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE TIMINT (J1,DT,EPS,WIL,NLEV,FIELD,FDT)
C--
C--   Aux. subr. TIMINT (J1,DT,EPS,WIL,NLEV,FIELD,FDT)
C--   Purpose : Perform time integration of FIELD at NLEV levels
C--             using tendency FDT
C--
      include "atparam.h"

      COMPLEX FIELD(MXNX,NLEV,2), FDT(MXNX,NLEV), FNEW(MXNX)

      EPS2=1.-2.*EPS

      IF (IX.EQ.IY*4) THEN
        DO K=1,NLEV
          CALL TRUNCT (FDT(1,K))
        ENDDO
      ENDIF

C the actual leap frog with the robert filter
      DO K=1,NLEV
        DO M=1,MXNX
          FNEW (M)     = FIELD(M,K,1)+DT*FDT(M,K)

          FIELD(M,K,1) = FIELD(M,K,J1) +  WIL*EPS*(FIELD(M,K,1)-
     &			 2*FIELD(M,K,J1)+FNEW(M))

C and here comes Williams' innovation to the filter
          FIELD(M,K,2) = FNEW(M)-(1-WIL)*EPS*(FIELD(M,K,1)-
     &			 2*FIELD(M,K,J1)+FNEW(M))

        ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE CGRATE (VOR,DIV,VORDT,DIVDT)
C--
C--   SUBROUTINE CGRATE (VOR,DIV,VORDT,DIVDT)
C--
C--   Purpose: Check growth rate of eddy kin. energy 
C--   Input  : VOR    = vorticity
C--            DIV    = divergence
C--            VORDT  = time derivative of VOR
C--            DIVDT  = time derivative of DIV
C--
      include "atparam.h"
      include "atparam1.h"

      COMPLEX VOR(MX,NX,KX), VORDT(MX,NX,KX), 
     &        DIV(MX,NX,KX), DIVDT(MX,NX,KX), TEMP(MX,NX)

      GRMAX=0.2/(86400.*2.)

      CDAMP=0.

      DO K=2,KX

        GRATE=0.
        RNORM=0.

        CALL INVLAP (VOR(1,1,K),TEMP)

        DO N=1,NX
         DO M=2,MX
           GRATE=GRATE-REAL(VORDT(M,N,K)*CONJG(TEMP(M,N)))
           RNORM=RNORM-REAL(  VOR(M,N,K)*CONJG(TEMP(M,N)))
         ENDDO
        ENDDO

        IF (GRATE.GT.GRMAX*RNORM) CDAMP =
     &      MAX(CDAMP,0.8*GRATE/RNORM)
C    &      MAX(CDAMP,(GRATE*GRATE)/(GRMAX*RNORM*RNORM))

      ENDDO

      IF (CDAMP.GT.0.) THEN

        print *, ' rot. wind damping enabled'

        DO K=1,KX
          DO N=1,NX
           DO M=2,MX
             VORDT(M,N,K)=VORDT(M,N,K)-CDAMP*VOR(M,N,K)
           ENDDO
          ENDDO
        ENDDO

      ENDIF


      CDAMP=0.

      DO K=2,KX

        GRATE=0.
        RNORM=0.

        CALL INVLAP (DIV(1,1,K),TEMP)

        DO N=1,NX
         DO M=2,MX
           GRATE=GRATE-REAL(DIVDT(M,N,K)*CONJG(TEMP(M,N)))
           RNORM=RNORM-REAL(  DIV(M,N,K)*CONJG(TEMP(M,N)))
         ENDDO
        ENDDO

        IF (GRATE.GT.GRMAX*RNORM) CDAMP =
     &      MAX(CDAMP,0.8*GRATE/RNORM)
C    &      MAX(CDAMP,(GRATE*GRATE)/(GRMAX*RNORM*RNORM))

      ENDDO

      IF (CDAMP.GT.0.) THEN

        print *, ' div. wind damping enabled'

        DO K=1,KX
          DO N=1,NX
           DO M=2,MX
             DIVDT(M,N,K)=DIVDT(M,N,K)-CDAMP*DIV(M,N,K)
           ENDDO
          ENDDO
        ENDDO

      ENDIF

C--
      RETURN
      END
