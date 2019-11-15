
      SUBROUTINE GRTEND (VORDT,DIVDT,TDT,PSDT,TRDT,J1,J2)
C--
C--   SUBROUTINE GRTEND (VORDT,DIVDT,TDT,PSDT,TRDT,J1,J2)
C--
C--   Purpose: compute non-linear tendencies in grid-point space
C              from dynamics and physical parametrizations,
C--            and convert them to spectral tendencies
C--
C--   dF/dt = T_dyn(F(J2)) + T_phy(F(J1))
C--
C--   Input:  J1 = time level index for physical tendencies 
C--           J2 = time level index for dynamical tendencies 
C--   Output: VORDT = spectral tendency of vorticity
C--           DIVDT = spectral tendency of divergence
C--           TDT   = spectral tendency of temperature
C--           PSDT  = spectral tendency of log(p_s)
C--           TRDT  = spectral tendency of tracers
C--

      include "atparam.h"
      include "atparam1.h"

      include "com_dyncon1.h"
      include "com_dyncon2.h"
      include "com_dynvar.h"

*** notes ****
c -- TG does not have to be computed at both time levels every time step,
c     I have left it this way to retain parallel structure with subroutine
c     using latitude loop
c -- memory can be reduced considerably eliminating TGG, computing VORG
c     only when needed, etc -- I have not optimized this subroutine for
c     routine use on the YMP
c -- results from grtend1.F should duplicate results from grtend.F
c                              -- Isaac
*************

      COMPLEX VORDT(MX,NX,KX), DIVDT(MX,NX,KX), TDT(MX,NX,KX),
     *        PSDT(MX,NX), TRDT(MX,NX,KX,NTR)

      COMPLEX DUMC(MX,NX,3), ZERO

      REAL UTEND(IX,IL,KX), VTEND(IX,IL,KX),
     *     TTEND(IX,IL,KX), TRTEND(IX,IL,KX,NTR)

      REAL UG(IX,IL,KX), VG(IX,IL,KX), TG(IX,IL,KX),
     *     TRG(IX,IL,KX,NTR),
     *     PX(IX,IL), PY(IX,IL),
     *     VORG(IX,IL,KX), DIVG(IX,IL,KX),
     *     SIGDT(IX,IL,KXP),
     *     PSTAR(IX,IL),
     *     UMEAN(IX,IL), VMEAN(IX,IL), DMEAN(IX,IL)

      REAL TGG(IX,IL,KX), PUV(IX,IL,KX),
     *     TEMP(IX,IL,KXP), SIGM(IX,IL,KXP), DUMR(IX,IL,3)


      ZERO=(0.,0.)

      iitest=0
      if (iitest.eq.1) print*,'inside GRTEND'

c -------------
c grid converts 

      if (iitest.eq.1) print*,'a'

      DO 10 K=1,KX

        CALL GRID(VOR(1,1,K,J2),VORG(1,1,K),1)
        CALL GRID(DIV(1,1,K,J2),DIVG(1,1,K),1)
        CALL GRID(  T(1,1,K,J2),  TG(1,1,K),1)

        DO 710 ITR=1,NTR
          CALL GRID(TR(1,1,K,J2,ITR),TRG(1,1,K,ITR),1)
  710   CONTINUE

        CALL UVSPEC(VOR(1,1,K,J2),DIV(1,1,K,J2),
     *              DUMC(1,1,1),DUMC(1,1,2))
        CALL GRID(DUMC(1,1,2),VG(1,1,K),2)
        CALL GRID(DUMC(1,1,1),UG(1,1,K),2)

        DO 11 J=1,IL
        DO 11 I=1,IX
          VORG(I,J,K)=VORG(I,J,K)+CORIOL(J)
   11   CONTINUE

   10 CONTINUE

      if (iitest.eq.1) print*,'b'
      DO 12 J=1,IL
      DO 12 I=1,IX
        UMEAN(I,J)=0.
        VMEAN(I,J)=0.
        DMEAN(I,J)=0.
   12 CONTINUE

      if (iitest.eq.1) print*,'c'
      DO 13 K=1,KX
        DO 113 J=1,IL
        DO 113 I=1,IX
          UMEAN(I,J)=UMEAN(I,J)+UG(I,J,K)*DHS(K)
          VMEAN(I,J)=VMEAN(I,J)+VG(I,J,K)*DHS(K)
          DMEAN(I,J)=DMEAN(I,J)+DIVG(I,J,K)*DHS(K)
  113   CONTINUE
   13 CONTINUE

c  compute tendency of log(surface pressure)

      if (iitest.eq.1) print*,'d'
C     PS(1,1,J2)=ZERO
      CALL GRAD(PS(1,1,J2),DUMC(1,1,2),DUMC(1,1,3))
      CALL GRID(DUMC(1,1,2),PX,2)
      CALL GRID(DUMC(1,1,3),PY,2)

      DO 14 J=1,IL
      DO 14 I=1,IX
        DUMR(I,J,1)=-UMEAN(I,J)*PX(I,J)-VMEAN(I,J)*PY(I,J)
   14 CONTINUE
      CALL SPEC(DUMR(1,1,1),PSDT)
      PSDT(1,1)=ZERO

c   compute "vertical" velocity  

      DO 16 J=1,IL
      DO 16 I=1,IX
        SIGDT(I,J,1)=0.
        SIGDT(I,J,KXP)=0.
        SIGM(I,J,1)=0.
        SIGM(I,J,KXP)=0.
   16 CONTINUE

c   (the following combination of terms is utilized later in the 
c       temperature equation)

      DO 250 K=1,KX
        DO 251 J=1,IL
        DO 251 I=1,IX
          PUV(I,J,K)=(UG(I,J,K)-UMEAN(I,J))*PX(I,J)+
     *     (VG(I,J,K)-VMEAN(I,J))*PY(I,J)
  251   CONTINUE
  250 CONTINUE
      if (iitest.eq.1) print*,'e'

      DO 17 K=1,KX
        DO 117 J=1,IL
        DO 117 I=1,IX
cspj    sigdt is the vertical velocity (in sigma coords)
          SIGDT(I,J,K+1)=SIGDT(I,J,K)
     *     -DHS(K)*(PUV(I,J,K)+DIVG(I,J,K)-DMEAN(I,J))
          SIGM(I,J,K+1)=SIGM(I,J,K)-DHS(K)*PUV(I,J,K)
  117   CONTINUE
   17 CONTINUE
 
c   subtract part of temperature field that is used as reference for 
c   implicit terms
      if (iitest.eq.1) print*,'f'

      DO 119 K=1,KX
      DO 118 J=1,IL
      DO 118 I=1,IX
        TGG(I,J,K)=TG(I,J,K)-TREF(K)
  118 CONTINUE
  119 CONTINUE

      DO 450 J=1,IL
      DO 450 I=1,IX
        PX(I,J)=RGAS*PX(I,J)
        PY(I,J)=RGAS*PY(I,J)
  450 CONTINUE

c   zonal wind tendency

      DO 399 J=1,IL
      DO 399 I=1,IX
        TEMP(I,J,1)=0.
        TEMP(I,J,KXP)=0.
  399 CONTINUE

      DO 400 K=2,KX
      DO 400 J=1,IL
      DO 400 I=1,IX
        TEMP(I,J,K)=SIGDT(I,J,K)*(UG(I,J,K)-UG(I,J,K-1))
  400 CONTINUE

      DO 500 K=1,KX
      DO 500 J=1,IL
      DO 500 I=1,IX
        UTEND(I,J,K)=VG(I,J,K)*VORG(I,J,K)
     *               -TGG(I,J,K)*PX(I,J)
     *               -(TEMP(I,J,K+1)+TEMP(I,J,K))*DHSR(K)
  500 CONTINUE

c  meridional wind tendency
      if (iitest.eq.1) print*,'g'

      DO 401 K=2,KX
      DO 401 J=1,IL
      DO 401 I=1,IX
        TEMP(I,J,K)=SIGDT(I,J,K)*(VG(I,J,K)-VG(I,J,K-1))
  401 CONTINUE

      DO 501 K=1,KX
      DO 501 J=1,IL
      DO 501 I=1,IX
        VTEND(I,J,K)=-UG(I,J,K)*VORG(I,J,K)
     *               -TGG(I,J,K)*PY(I,J)
     *               -(TEMP(I,J,K+1)+TEMP(I,J,K))*DHSR(K)
  501 CONTINUE

c  temperature tendency

      DO 402 K=2,KX
      DO 402 J=1,IL
      DO 402 I=1,IX
        TEMP(I,J,K)=SIGDT(I,J,K)*(TGG(I,J,K)-TGG(I,J,K-1))
     *             +SIGM(I,J,K)*(TREF(K)-TREF(K-1))
  402 CONTINUE

      DO 502 K=1,KX
      DO 502 J=1,IL
      DO 502 I=1,IX
        TTEND(I,J,K)=TGG(I,J,K)*DIVG(I,J,K)
     *               -(TEMP(I,J,K+1)+TEMP(I,J,K))*DHSR(K)
     *               +FSGR(K)*TGG(I,J,K)*(SIGDT(I,J,K+1)+SIGDT(I,J,K))
     *               +TREF3(K)*(SIGM(I,J,K+1)+SIGM(I,J,K))
     *               +AKAP*(TG(I,J,K)*PUV(I,J,K)
     *               -TGG(I,J,K)*DMEAN(I,J))
  502 CONTINUE

      if (iitest.eq.1) print*,'h'
c  tracer tendency

      DO 443 ITR=1,NTR 

        DO 522 K=2,KX
        DO 522 J=1,IL
        DO 522 I=1,IX
          TEMP(I,J,K)=SIGDT(I,J,K)*(TRG(I,J,K,ITR)-TRG(I,J,K-1,ITR)) 
  522 CONTINUE


cspj for moisture, vertical advection is not possible between top 
cspj two layers
ckuch three layers
c      if(iinewtrace.eq.1)then
        DO K=2,3
         DO J=1,IL
          DO I=1,IX
           TEMP(I,J,K)=0.
          ENDDO
         ENDDO
        ENDDO
c      endif

        DO 542 K=1,KX
        DO 542 J=1,IL
        DO 542 I=1,IX
          TRTEND(I,J,K,ITR)=TRG(I,J,K,ITR)*DIVG(I,J,K)
     *                      -(TEMP(I,J,K+1)+TEMP(I,J,K))*DHSR(K)
  542 CONTINUE

  443 CONTINUE

      if (iitest.eq.1) print*,'h'

******************** physics ****************************

      CALL GEOP (J1)

      CALL PHYPAR (VOR(1,1,1,J1),DIV(1,1,1,J1),T(1,1,1,J1),
     *             TR(1,1,1,J1,1),PHI,PS(1,1,J1),
     *             UTEND,VTEND,TTEND,TRTEND)


*********************************************************
      if (iitest.eq.1) print*,'i'

      DO 18 K=1,KX

c  convert u and v tendencies to vor and div spectral tendencies
c  vdspec takes a grid u and a grid v and converts them to 
c  spectral vor and div

        CALL VDSPEC(UTEND(1,1,K),VTEND(1,1,K),VORDT(1,1,K),
     *              DIVDT(1,1,K),2)

c  add lapl(0.5*(u**2+v**2)) to div tendency,
c  and add div(vT) to spectral t tendency

        DO 20 j=1,il
        DO 20 I=1,IX
          DUMR(I,J,1)=0.5*(UG(I,J,K)*UG(I,J,K)
     *                    +VG(I,J,K)*VG(I,J,K))
          DUMR(I,J,2)=-UG(I,J,K)*TGG(I,J,K)
          DUMR(I,J,3)=-VG(I,J,K)*TGG(I,J,K)
   20   CONTINUE

c  divergence tendency

        CALL SPEC(DUMR(1,1,1),DUMC(1,1,1))
        CALL LAP (DUMC(1,1,1),DUMC(1,1,2))

cfk--   Change to keep dimensions 
      
        DO 222 M=1,MX
        DO 222 N=1,NX
          DIVDT(M,N,K)=DIVDT(M,N,K)-DUMC(M,N,2)
  222   CONTINUE

c  temperature tendency

        CALL VDSPEC(DUMR(1,1,2),DUMR(1,1,3),DUMC(1,1,1),
     *              TDT(1,1,K),2)
        CALL SPEC(TTEND(1,1,K),DUMC(1,1,2))

cfk--   Change to keep dimensions 

        DO 22 M=1,MX
        DO 22 N=1,NX
          TDT(M,N,K)=TDT(M,N,K)+DUMC(M,N,2)
   22   CONTINUE

c tracer tendency

        DO 31 ITR=1,NTR

          DO 30 J=1,IL
          DO 30 I=1,IX
            DUMR(I,J,2)=-UG(I,J,K)*TRG(I,J,K,ITR)
            DUMR(I,J,3)=-VG(I,J,K)*TRG(I,J,K,ITR)
   30     CONTINUE
c 
          CALL SPEC(TRTEND(1,1,K,ITR),DUMC(1,1,2))
          CALL VDSPEC(DUMR(1,1,2),DUMR(1,1,3),DUMC(1,1,1),
     *                TRDT(1,1,K,ITR),2)

cfk--   Change to keep dimensions 

          DO 32 M=1,MX
          DO 32 N=1,NX
            TRDT(M,N,K,ITR)=TRDT(M,N,K,ITR)+DUMC(M,N,2)
   32     CONTINUE

   31   CONTINUE

   18 CONTINUE
      
      RETURN
      END 
