
      SUBROUTINE FORINT (NGP,IMON,FMON,FORMNT,FOR1,MONTHS)  

C--   Aux. routine FORINT : linear interpolation of monthly-mean forcing
C--IO sx 12 months in a forcing cycle

      REAL FORMNT(NGP,*), FOR1(NGP)

      IF (FMON.LE.0.5) THEN
        IMON2 = IMON-1
        IF (IMON.EQ.1) IMON2 = MONTHS
        WMON = 0.5-FMON
      ELSE
        IMON2 = IMON+1
        IF (IMON.EQ.MONTHS) IMON2 = 1
        WMON = FMON-0.5
      ENDIF

      DO J=1,NGP
        FOR1(J) = FORMNT(J,IMON)+WMON*(FORMNT(J,IMON2)-FORMNT(J,IMON))
      ENDDO
C--
      RETURN
      END

      subroutine FORIN5 (ngp,imon,fmon,formnt,for1,months)

C--   Aux. routine FORIN5 : non-linear, mean-conserving interpolation 
C--                         of monthly-mean forcing fields
C--IO sx 12 months in a forcing cycle

      real formnt(ngp,months), for1(ngp)

      im2 = imon-2
      im1 = imon-1
      ip1 = imon+1
      ip2 = imon+2

      if (im2.lt.1)  im2 = im2+months
      if (im1.lt.1)  im1 = im1+months
      if (ip1.gt.months) ip1 = ip1-months
      if (ip2.gt.months) ip2 = ip2-months
 
      c0 = 1./real(months)
      t0 = c0*fmon
      t1 = c0*(1.-fmon)
      t2 = 0.25*fmon*(1-fmon)

      wm2 =        -t1   +t2
      wm1 =  -c0 +8*t1 -6*t2
      w0  = 7*c0      +10*t2     
      wp1 =  -c0 +8*t0 -6*t2
      wp2 =        -t0   +t2 

      do j=1,ngp
        for1(j) = wm2*formnt(j,im2)+wm1*formnt(j,im1)
     &           + w0*formnt(j,imon)
     &           +wp1*formnt(j,ip1)+wp2*formnt(j,ip2)
      enddo

      return
      end

