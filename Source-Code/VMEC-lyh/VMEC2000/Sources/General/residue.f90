      SUBROUTINE residue_par (gcr, gcz, gcl)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: rss, zcs, rsc, zcc,                               &
                             meven, modd, ntmax, signgs
      USE realspace, ONLY: phip
      USE xstuff
      USE precon2d
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: tlglob_arr, trglob_arr,                  &
                                      lactive, SAXLASTNTYPE
      USE blocktridiagonalsolver, ONLY: L_COLSCALE

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntmax), INTENT(INOUT) ::           &
         gcr, gcz, gcl
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n0=0, m0=0, m1=1
      INTEGER, PARAMETER :: n3d=0, nasym=1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nsfix, jedge, delIter
      REAL(dp) :: r1, fac, totals(3), rznorm, edge(1)
      INTEGER :: l
!-----------------------------------------------
      CALL second0 (treson)

#ifdef _HBANGLE
!FREE-BDY RFP MAY NEED THIS TO IMPROVE CONVERGENCE (SPH 022514)
      IF (lfreeb .AND. lrfp) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) THEN
            fac = 1.E-1_dp
         END IF
         gcz(ns,0,m0,:) = fac*gcz(ns,0,m0,:)
      END IF
#else
!
!     SYMMETRIC PERTURBATIONS (BASED ON POLAR RELATIONS):
!        Rss(n) = Zcs(n), n != 0
!     ASYMMETRIC PERTURBATIONS:
!        Rsc(n) = Zcc(n), ALL n
!
!     INTERNALLY:
!        XC(rss) = .5*(Rss + Zcs), XC(zcs) = .5*(Rss - Zcs) -> 0
!        XC(rsc) = .5*(Rsc + Zcc), XC(zcc) = .5*(Rsc - Zcc) -> 0
!     THIS IMPLIES THE CONSTRAINT
!        3D ONLY : GC(zcs) = 0;  
!        ASYM:     GC(zcc) = 0
!

      IF (lthreed) THEN
         CALL constrain_m1_par(gcr(:,m1,:,rss), gcz(:,m1,:,zcs))
      END IF
      IF (lasym) THEN
         CALL constrain_m1_par(gcr(:,m1,:,rsc), gcz(:,m1,:,zcc))
      END IF

      IF (lfreeb .AND. lrfp) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) THEN
            fac = 1.E-1_dp
         END IF
         gcr(0,m0,ns,:) = fac*gcr(0,m0,ns,:)
         gcz(0,m0,ns,:) = fac*gcz(0,m0,ns,:)
      END IF
#endif

!     PRECONDITIONER MUST BE CALCULATED USING RAW (UNPRECONDITIONED) FORCES
      IF (ictrl_prec2d .GE. 2 .OR. ictrl_prec2d .EQ. -1) THEN
         RETURN
      END IF

!
!     COMPUTE INVARIANT RESIDUALS
!
      r1 = one/(2*r0scale)**2
      jedge = 0
      delIter = iter2-iter1

      IF (delIter       .lt. 50 .and.                                        &
          (fsqr + fsqz) .LT. 1.E-6_dp) THEN
         jedge = 1
      ENDIF

      rznorm = r1*fnorm
      CALL getfsq_local(gcr, gcz, totals(1:2), jedge)
      totals(3) = SUM(gcl(:,:,tlglob:trglob,:)**2)
      IF(rank .EQ. nranks-1) THEN
         fedge = rznorm*SUM(gcr(:,:,ns,:)**2 + gcz(:,:,ns,:)**2)
      END IF
      CALL sum_surface_scalars(totals, 3)
      fsqr = rznorm*totals(1)
      fsqz = rznorm*totals(2)
      fsql = fnormL*totals(3)
      IF (l_coupled_gmres) THEN
         gmres_merit = fsqr + fsqz + fsql
         IF (lfreeb .AND. jedge.EQ.0) THEN
            edge = 0
            IF (rank.EQ.nranks-1) edge = fedge
            CALL sum_surface_scalars(edge, 1)
            gmres_merit = gmres_merit + edge(1)
         END IF
      END IF
!
!     PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!
      IF (ictrl_prec2d .EQ. 1) THEN
         
         IF (l_colscale .AND. lactive) THEN
            CALL SAXLASTNTYPE(pgc, pcol_scale, pgc)
         END IF

         LRESIDUECALL = .TRUE.
         CALL block_precond_par(pgc)
         LRESIDUECALL = .FALSE.

         IF (.NOT.lfreeb .AND. ANY(gcr(:,:,ns,:) .NE. zero)) THEN
            STOP 'gcr(ns) != 0 for fixed boundary in residue'
         END IF
         IF (.NOT.lfreeb .AND. ANY(gcz(:,:,ns,:) .NE. zero)) THEN
            STOP 'gcz(ns) != 0 for fixed boundary in residue'
         END IF
         IF (ANY(gcl(1:,m0,:,zsc) .NE. zero)) THEN
            STOP 'gcl(m=0,n>0,sc) != 0 in residue'
         END IF
         IF (lthreed) THEN
            IF (ANY(gcl(n0,:,:,zcs).NE.zero)) THEN
               STOP 'gcl(n=0,m,cs) != 0 in residue'
            END IF
         END IF

         fsqr1 = SUM(gcr*gcr)
         fsqz1 = SUM(gcz*gcz)
         fsql1 = SUM(gcl*gcl)

      ELSE
!        m = 1 constraint scaling

         IF (lthreed) THEN
            CALL scale_m1_par(gcr(:,m1,:,rss), gcz(:,m1,:,zcs))
         END IF
         IF (lasym) THEN
            CALL scale_m1_par(gcr(:,m1,:,rsc), gcz(:,m1,:,zcc))
         END IF

         jedge = 0
         CALL scalfor_par (gcr, arm, brm, ard, brd, crd, jedge)
         jedge = 1
         CALL scalfor_par (gcz, azm, bzm, azd, bzd, crd, jedge)

         CALL getfsq_local(gcr, gcz, totals(1:2), m1)
         DO l = tlglob, trglob
            gcl(:,:,l,:) = pfaclam(:,:,l,:)*gcl(:,:,l,:)
         END DO
         totals(3) = 0
         DO l = MAX(2,tlglob), trglob
            totals(3) = totals(3) + SUM(gcl(:,:,l,:)**2)
         END DO
         CALL sum_surface_scalars(totals, 3)
         fsqr1 = fnorm1*totals(1)
         fsqz1 = fnorm1*totals(2)
         fsql1 = hs*totals(3)

         CALL PadSides(pgc)  

      ENDIF

      CALL second0 (tresoff)
      residue_time = residue_time + (tresoff-treson)

      END SUBROUTINE residue_par

      SUBROUTINE constrain_m1_par(gcr, gcz)
      USE vmec_main
      USE parallel_include_module
      USE precon2d, ONLY: m1_force_constrained
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(INOUT) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: temp
!-----------------------------------------------
!
!     COMPUTE INTERNAL gr, gz
!     NOTE: internal gz => 0 for both values of lconm1 (although gz is different)
!     FOR lconm1=T, gcr(internal) = gcr+gcz, gcz(internal) = gcr-gcz->0
!
      ALLOCATE(temp(0:ntor,ns))
      IF (lconm1) THEN
         temp(:,tlglob:trglob) = gcr(:,tlglob:trglob)
         gcr(:,tlglob:trglob) = osqrt2*(gcr(:,tlglob:trglob) +                 &
                                        gcz(:,tlglob:trglob))
         gcz(:,tlglob:trglob) = osqrt2*(temp(:,tlglob:trglob) -                &
                                        gcz(:,tlglob:trglob))
      END IF

!v8.50: ADD iter2<2 so reset=<WOUT_FILE> works
      IF (m1_force_constrained()) THEN
         gcz(:,tlglob:trglob) = 0
      END IF
 
      DEALLOCATE(temp)
      END SUBROUTINE constrain_m1_par

      SUBROUTINE scale_m1_par(gcr, gcz)
      USE vmec_main
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,ns), INTENT(inout) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: nodd=2
      INTEGER :: n
      REAL(dp) :: fac(ns)
!-----------------------------------------------
      IF (.not.lconm1) RETURN

      fac(tlglob:trglob) = (ard(tlglob:trglob,nodd) +                          &
                            brd(tlglob:trglob,nodd))                           &
                         / (ard(tlglob:trglob,nodd) +                          &
                            brd(tlglob:trglob,nodd) +                          &
                            azd(tlglob:trglob,nodd) +                          &
                            bzd(tlglob:trglob,nodd))
      DO n = 0, ntor
         gcr(n,tlglob:trglob) = fac(tlglob:trglob)*gcr(n,tlglob:trglob)
      END DO

      fac(tlglob:trglob) = (azd(tlglob:trglob,nodd) +                          &
                            bzd(tlglob:trglob,nodd))                           &
                         / (ard(tlglob:trglob,nodd) +                          &
                            brd(tlglob:trglob,nodd) +                          &
                            azd(tlglob:trglob,nodd) +                          &
                            bzd(tlglob:trglob,nodd))
      DO n = 0, ntor
         gcz(n,tlglob:trglob) = fac(tlglob:trglob)*gcz(n,tlglob:trglob)
      END DO
 
      END SUBROUTINE scale_m1_par

      SUBROUTINE residue (gcr, gcz, gcl)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: rss, zcs, rsc, zcc,                               &
                             meven, modd, ntmax, signgs
      USE realspace, ONLY: phip
      USE xstuff
      USE precon2d
#ifdef _HBANGLE
      USE angle_constraints, ONLY: scalfor_rho
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: &
        gcr, gcz, gcl
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n0=0, m0=0, m1=1
      INTEGER, PARAMETER :: n3d=0, nasym=1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nsfix, jedge, delIter
      REAL(dp) :: r1, fac
!-----------------------------------------------
!
!     IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!     INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!     (ZCS = RSS, ZSS = RCS ARE THE CORRECT POLAR RELATIONS)
!

#ifdef _HBANGLE
!FREE-BDY RFP MAY NEED THIS TO IMPROVE CONVERGENCE (SPH 022514)
      IF (lfreeb .AND. lrfp) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) THEN
            fac = 1.E-1_dp
         END IF
         gcz(ns,0,m0,:) = fac*gcz(ns,0,m0,:)
      END IF
#else
!
!     SYMMETRIC PERTURBATIONS (BASED ON POLAR RELATIONS):
!        RSS(n) = ZCS(n), n != 0
!     ASYMMETRIC PERTURBATIONS:
!        RSC(n) = ZCC(n), ALL n
!
!     INTERNALLY:
!        XC(rss) = .5*(Rss + Zcs), XC(zcs) = .5*(Rss - Zcs) -> 0
!        XC(rsc) = .5*(Rsc + Zcc), XC(zcc) = .5*(Rsc - Zcc) -> 0
!     THIS IMPLIES THE CONSTRAINT
!        3D ONLY : GC(zcs) = 0;  
!        ASYM:     GC(zcc) = 0
!

      IF (lthreed) THEN
         CALL constrain_m1(gcr(:,:,m1,rss), gcz(:,:,m1,zcs))
      END IF
      IF (lasym) THEN
         CALL constrain_m1(gcr(:,:,m1,rsc), gcz(:,:,m1,zcc))
      END IF

!FREE-BDY RFP MAY NEED THIS TO IMPROVE CONVERGENCE (SPH 022514)
      IF (lfreeb .AND. lrfp) THEN
         fac = 0
         IF (ictrl_prec2d .EQ. 0) THEN
            fac = 1.E-1_dp
         END IF
         gcr(ns,0,m0,:) = fac*gcr(ns,0,m0,:)
         gcz(ns,0,m0,:) = fac*gcz(ns,0,m0,:)
      END IF
#endif

!     PRECONDITIONER MUST BE CALCULATED USING RAW (UNPRECONDITIONED) FORCES
      IF (ictrl_prec2d .GE. 2 .OR.                                             &
          ictrl_prec2d .EQ. -1) THEN
         RETURN
      END IF

!
!     COMPUTE INVARIANT RESIDUALS
!
      r1 = one/(2*r0scale)**2
      jedge = 0
!SPH-JAH013108: MUST INCLUDE EDGE FORCE (INITIALLY) FOR V3FITA TO WORK
!ADD A V3FIT RELATED FLAG? ADD fsq criterion first
      delIter = iter2 - iter1

      IF (delIter       .lt. 50 .and.                                          &
          (fsqr + fsqz) .lt. 1.E-6_dp) THEN
         jedge = 1
      END IF

      CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

      fsql = fnormL*SUM(gcl*gcl)
      fedge = r1*fnorm*SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)
      IF (l_coupled_gmres) THEN
         gmres_merit = fsqr + fsqz + fsql
         IF (lfreeb .AND. jedge.EQ.0) gmres_merit = gmres_merit + fedge
      END IF

!
!     PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!

      IF (ictrl_prec2d .EQ. 1) THEN

         CALL block_precond(gc)

         IF (.not.lfreeb .and. ANY(gcr(ns,:,:,:) .ne. zero)) THEN
            STOP 'gcr(ns) != 0 for fixed boundary in residue'
         END IF
         IF (.not.lfreeb .and. ANY(gcz(ns,:,:,:) .ne. zero)) THEN
            STOP 'gcz(ns) != 0 for fixed boundary in residue'
         END IF
         IF (ANY(gcl(:,1:,0,zsc) .ne. zero)) THEN
            STOP 'gcl(m=0,n>0,sc) != 0 in residue'
         END IF
         IF (lthreed) THEN
            IF (ANY(gcl(:,n0,:,zcs).NE.zero)) THEN
               STOP 'gcl(n=0,m,cs) != 0 in residue'
            END IF
         END IF

         fsqr1 = SUM(gcr*gcr)
         fsqz1 = SUM(gcz*gcz)
         fsql1 = SUM(gcl*gcl)

      ELSE
#ifdef _HBANGLE
         CALL scalfor_rho(gcr, gcz)
#else
!        m = 1 constraint scaling
         IF (lthreed) THEN
            CALL scale_m1(gcr(:,:,m1,rss), gcz(:,:,m1,zcs))
         END IF
         IF (lasym) THEN
            CALL scale_m1(gcr(:,:,m1,rsc), gcz(:,:,m1,zcc))
         END IF

         jedge = 0
         CALL scalfor (gcr, arm, brm, ard, brd, crd, jedge)
         jedge = 1
         CALL scalfor (gcz, azm, bzm, azd, bzd, crd, jedge)
#endif

!SPH: add fnorm1 ~ 1/R**2, since preconditioned forces gcr,gcz ~ Rmn or Zmn
         CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, m1)
#ifdef _HBANGLE
!
!     TO IMPROVE CONVERGENCE, REDUCE FORCES INITIALLY IF THEY ARE TOO LARGE
!
         fac = .5_dp
         IF ((iter2 - iter1) .LT. 25 .AND.                                     &
             (fsqr + fsqz)   .GT. 1.E-2_dp) THEN
            fac = fac / SQRT(1.E2_dp*(fsqr+fsqz))
         END IF
         gcr = fac*gcr
         gcz = fac*gcz
#endif

!SPH: THIS IS NOT INVARIANT UNDER PHIP->A*PHIP, AM->A**2*AM IN PROFIL1D
!     (EXTCUR -> A*EXTCUR for FREE BOUNDARY)

         gcl = faclam*gcl
         fsql1 = hs*SUM(gcl*gcl)
!030514      fsql1 = hs*lamscale**2*SUM(gcl*gcl)

      END IF

      END SUBROUTINE residue

      SUBROUTINE constrain_m1(gcr, gcz)
      USE vmec_main
      USE precon2d, ONLY: m1_force_constrained
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor), INTENT(INOUT) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp) :: temp(ns,0:ntor)
!-----------------------------------------------
!
!     COMPUTE INTERNAL gr, gz
!     NOTE: internal gz => 0 for both values of lconm1 (although gz is different)
!     FOR lconm1=T, gcr(internal) = gcr+gcz, gcz(internal) = gcr-gcz->0
!
      IF (lconm1) THEN
         temp = gcr
         gcr = osqrt2*(gcr + gcz)
         gcz = osqrt2*(temp - gcz)
      END IF

!v8.50: ADD iter2<2 so reset=<WOUT_FILE> works
      IF (m1_force_constrained()) THEN
         gcz = 0
      END IF
 
      END SUBROUTINE constrain_m1

      SUBROUTINE scale_m1(gcr, gcz)
      USE vmec_main
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: nodd=2
      INTEGER :: n
      REAL(dp) :: fac(ns)
!-----------------------------------------------
      IF (.not.lconm1) RETURN

      fac = (ard(1:ns,nodd) + brd(1:ns,nodd))                                        &
          / (ard(1:ns,nodd) + brd(1:ns,nodd) + azd(1:ns,nodd) + bzd(1:ns,nodd))
      DO n = 0, ntor
         gcr(:,n) = fac*gcr(:,n)
      END DO

      fac = (azd(1:ns,nodd) + bzd(1:ns,nodd))                                        &
          / (ard(1:ns,nodd) + brd(1:ns,nodd) + azd(1:ns,nodd) + bzd(1:ns,nodd))
      DO n = 0, ntor
         gcz(:,n) = fac*gcz(:,n)
      END DO
 
      END SUBROUTINE scale_m1

