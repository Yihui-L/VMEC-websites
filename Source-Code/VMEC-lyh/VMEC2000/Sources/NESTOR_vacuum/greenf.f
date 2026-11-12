      RECURSIVE SUBROUTINE greenf (delgr, delgrp, ip, gsave, dsave)
      USE vacmod
      USE vparams, ONLY: one
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ip
      REAL(dp), DIMENSION(nuv), INTENT(OUT) :: delgr, delgrp
      REAL(dp), DIMENSION(nuv), INTENT(OUT) :: gsave, dsave
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, DIMENSION(2) :: ilow, ihigh
      INTEGER :: ivoff, iskip, iuoff, i, kp, nloop
      REAL(dp) :: ftemp, htemp, ga1, ga2
      REAL(dp) :: xip, yip, xper, yper, sxsave, sysave
C-----------------------------------------------
!
!     ON ENTRANCE, IP IS THE INDEX OF THE PRIMED MESH POINT (lies in 1st field period)
!
!     ON EXIT, DELGR IS THE DIFFERENCE OF "GREEN'S FUNCTION"
!     AND ANALYTIC APPROXIMATION, SUMMED OVER ALL FIELD PERIODS
!     DELGRP IS DIFFERENCE OF DERIVATIVE OF "GREEN'S FUNCTION"
!     AND ANALYTIC APPROXIMATION.
!
!     BOTH THESE QUANTITIES ARE COMPUTED FOR ALL UNPRIMED U,V POINTS IN ONE FIELD PERIOD,
!     FOR THIS FIXED PRIMED POINT (IP).
!
!
!     COMPUTE OFFSETS FOR U,V ANGLE DIFFERENCES AND CONSTANTS
!
      ilow(1) = 1
      ilow(2) = ip + 1
      ihigh(1) = ip - 1
      ihigh(2) = nuv
      ivoff = nuv + 1 - ip
      iskip = (ip - 1)/nv
      iuoff = nuv - nv*iskip
      xip = rcosuv(ip)             !x == r*COS(ip), in 1st field period
      yip = rsinuv(ip)             !y == r*SIN(ip), in 1st field period
      delgr  = 0
      delgrp = 0

!
!     COMPUTE FIELD-PERIOD INVARIANT VECTORS
!
!     NOTE: |x - x'|**2 = gsave - 2*[x*x' + y*y']
!
      DO i = 1, nuv
        gsave(i) = rzb2(ip) + rzb2(i) - 2*z1b(ip)*z1b(i)
        dsave(i) = drv(ip) + z1b(i)*snz(ip)
      END DO

!
!     SUM OVER FIELD-PERIODS (NVPER=NFPER) OR INTEGRATE OVER NV (NVPER=64) IF NV == 1
!
!     NOTE THE SURFACE NORMAL SNORM == Xu cross Xv = NP*[SNR, SNV, SNZ]
!     IS PERIODIC ON EACH FIELD PERIOD
!
      DO kp = 1, nvper
         xper = xip*cosper(kp) - yip*sinper(kp)              !x(ip) in field period kp
         yper = yip*cosper(kp) + xip*sinper(kp)              !y(ip) in field period kp
         sxsave = (snr(ip)*xper - snv(ip)*yper)/r1b(ip)
         sysave = (snr(ip)*yper + snv(ip)*xper)/r1b(ip)

         IF (kp.EQ.1 .OR. nv.EQ.1) THEN

!        Form and use the analytic terms at each nonsingular point.
            DO nloop = 1, 2
               IF (kp.GT.1 .AND. nloop.EQ.2) CYCLE
               DO i = ilow(nloop), ihigh(nloop)
                 ga1 = tanu(i+iuoff)*(guu_b(ip)*tanu(i+iuoff)
     &                              + guv_b(ip)*tanv(i+ivoff))
     &               + gvv_b(ip)*tanv(i+ivoff)*tanv(i+ivoff)
                 ga2 = tanu(i+iuoff)*(auu(ip)*tanu(i+iuoff)
     &                              + auv(ip)*tanv(i+ivoff))
     &               + avv(ip)*tanv(i+ivoff)*tanv(i+ivoff)
                 ga2 = ga2/ga1
                 ga1 = one/SQRT(ga1)
                 ftemp = one/(gsave(i)
     &                 - 2*(xper*rcosuv(i) + yper*rsinuv(i)))
                 htemp = SQRT(ftemp)
                 delgrp(i) = delgrp(i) - ga2*ga1 + ftemp*htemp*
     &              (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
                 delgr(i) = delgr(i) + htemp - ga1
               END DO
            END DO

            IF (kp.EQ.nvper .AND. nv.EQ.1) THEN
                delgrp = delgrp/nvper
                delgr  = delgr /nvper
            END IF

            ivoff = ivoff + 2*nu
            ihigh(1) = nuv

         ELSE
            DO i = 1, nuv
              ftemp = one/(gsave(i)
     1                 -  2*(xper*rcosuv(i) + yper*rsinuv(i)))
              htemp = SQRT(ftemp)
              delgrp(i) = delgrp(i) + ftemp*htemp*
     1           (rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
              delgr(i) = delgr(i) + htemp
           END DO
         ENDIF
      END DO

      END SUBROUTINE greenf

