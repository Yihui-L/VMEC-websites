      MODULE gmres_mod
      USE vmec_main, ONLY: dp, rprec, neqs, ns, nthreed, 
     1                     one, fsqr, fsqz, fsql
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: CopyLastNtype, SaxpbyLastNtype,
     1                                SaxpbyLastNs, Saxpby1LastNs,
     2                                GetDerivLastNs, SaxLastNtype
      USE precon2d, ONLY: ictrl_prec2d,
     &                   trial_context => l_gmres_trial, gmres_merit,
     &                   l_coupled_gmres
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
      IMPLICIT NONE
      INTEGER :: nfcn = 0, lqmrstep = 0
      INTEGER :: ier_flag_res
      INTEGER, PRIVATE :: legacy_matvec_status = 0
      LOGICAL :: lqmr, lfirst
      LOGICAL, PARAMETER :: lscreen0 = .FALSE.
!     A nonlinear trial must not advance the outer vacuum schedule or
!     change the history-dependent native residual mask between calls.
      PRIVATE :: trial_context
      LOGICAL, PRIVATE :: bad_matvec = .FALSE.
      INTEGER, PRIVATE :: saved_irst, saved_ivac, saved_nvacskip,
     &                    saved_nvskip0
      REAL(dp), PRIVATE :: saved_fsq(6)

!
!     nfcn :  number of calls to function (funct3d)
!     lqmr :  logical, used by external programs to control calling these routines
!
      CONTAINS

      SUBROUTINE begin_gmres_trials
      USE vmec_main, ONLY: irst, ivac, nvacskip, nvskip0,
     &                     fsqr1, fsqz1, fsql1
      saved_irst = irst
      saved_ivac = ivac
      saved_nvacskip = nvacskip
      saved_nvskip0 = nvskip0
      saved_fsq = (/fsqr, fsqz, fsql, fsqr1, fsqz1, fsql1/)
      trial_context = .TRUE.
      bad_matvec = .FALSE.
      END SUBROUTINE begin_gmres_trials

      SUBROUTINE prepare_gmres_trial
      USE vmec_main, ONLY: irst, ivac, nvacskip, nvskip0,
     &                     fsqr1, fsqz1, fsql1
      ier_flag_res = 0
      IF (.NOT.trial_context) RETURN
      irst = saved_irst
      ivac = saved_ivac
      nvacskip = saved_nvacskip
      nvskip0 = saved_nvskip0
      fsqr = saved_fsq(1); fsqz = saved_fsq(2)
      fsql = saved_fsq(3); fsqr1 = saved_fsq(4)
      fsqz1 = saved_fsq(5); fsql1 = saved_fsq(6)
      END SUBROUTINE prepare_gmres_trial

      SUBROUTINE end_gmres_trials
      USE vmec_main, ONLY: ivac, nvacskip, nvskip0
!     Keep the force/diagnostics of the accepted or restored physical
!     state. Only undo scheduling changes caused by inner evaluations.
      ivac = saved_ivac
      nvacskip = saved_nvacskip
      nvskip0 = saved_nvskip0
      trial_context = .FALSE.
      END SUBROUTINE end_gmres_trials

      SUBROUTINE collective_trial_ok(ok, parallel)
      LOGICAL, INTENT(INOUT) :: ok
      LOGICAL, INTENT(IN) :: parallel
#if defined(MPI_OPT)
      LOGICAL :: all_ok
      IF (parallel) THEN
         CALL MPI_ALLREDUCE(ok, all_ok, 1, MPI_LOGICAL, MPI_LAND,
     &                      RUNVMEC_COMM_WORLD, MPI_ERR)
         ok = all_ok
      END IF
#endif
      END SUBROUTINE collective_trial_ok

      LOGICAL FUNCTION valid_gmres_trial(parallel) RESULT(ok)
      USE vmec_main, ONLY: irst, fsqr1, fsqz1, fsql1
      LOGICAL, INTENT(IN) :: parallel
      REAL(dp) :: values(6)
      values = (/fsqr, fsqz, fsql, fsqr1, fsqz1, fsql1/)
      ok = ier_flag_res.EQ.0 .AND. irst.EQ.1
      ok = ok .AND. ALL(ieee_is_finite(values))
      IF (ok) ok = ALL(values.GE.0._dp)
      CALL collective_trial_ok(ok, parallel)
      END FUNCTION valid_gmres_trial

      SUBROUTINE matvec_par (ploc, Ap, nloc)
      USE blocktridiagonalsolver, ONLY: ParMatVec
      USE stel_kinds
      USE xstuff, ONLY: pxc, px0=>pxsave, pgc0=>pxcdot, pgc
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(IN)   :: nloc
      REAL(dp), INTENT(IN)  :: 
     &   ploc(ntmaxblocksize,tlglob:trglob)
      REAL(dp), INTENT(OUT) :: 
     &   Ap(ntmaxblocksize,tlglob:trglob)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: mblk_size, istat
      REAL(dp) :: delta, lmax, gmax, pmax, apmax
      LOGICAL :: ok
C-----------------------------------------------
      mblk_size = neqs/ns
!
!     Computes linearized matrix product A*p = [F(x0+delta*p) - F(x0)]/delta, about point x0
!     Must scale p so delta*max|p| ~ sqrt(epsilon) to get accurate numerical derivative
!     Because the block preconditioner is applied in funct3d, the net result
!     should be Ap ~ -p.

      delta = SQRT(EPSILON(delta))

      LACTIVE0: IF (lactive) THEN
         IF (nloc .NE. (trglob - tlglob + 1)*mblk_size) THEN
            STOP 'nloc wrong in matvec_par'
         END IF
         lmax = SUM(ploc*ploc)
         CALL MPI_ALLREDUCE(lmax, pmax, 1, MPI_REAL8, MPI_SUM, NS_COMM,
     &                      MPI_ERR)
         pmax = SQRT(pmax)
         delta = delta/MAX(delta, pmax)

         CALL SaxpbyLastNs(delta, ploc, one, px0, pxc)

         CALL last_ntype_par
         CALL PadSides(pxc)
      END IF LACTIVE0

      CALL prepare_gmres_trial
      CALL funct3d_par(lscreen0, ier_flag_res)
      
      IF (lactive) THEN
         CALL last_ns_par
         CALL GetDerivLastNs(pgc, pgc0, delta, Ap)
      ENDIF

      IF (trial_context) THEN
         ok = valid_gmres_trial(.TRUE.)
         IF (lactive) ok = ok .AND. ALL(ieee_is_finite(Ap))
         CALL collective_trial_ok(ok, .TRUE.)
         IF (.NOT.ok) THEN
!           Let the bounded inner solver return, then reject its step.
            bad_matvec = .TRUE.
            IF (lactive) Ap = -ploc
         END IF
      ELSE IF (ier_flag_res.NE.0 .AND. rank.EQ.0) THEN
         PRINT *,'IN MATVEC_PAR, IER_FLAG = ', ier_flag_res
      END IF

 90   CONTINUE
      nfcn = nfcn + 1

      END SUBROUTINE matvec_par

      SUBROUTINE GetNLForce_par(xcstate, fsq_nl, bnorm)
      USE xstuff, ONLY: pxc, pgc, x0=>pxsave
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: xcstate(neqs), bnorm
      REAL(dp), INTENT(OUT) :: fsq_nl
!-----------------------------------------------
!undo internal gmres normalization

      LACTIVE0: IF (lactive) THEN
         CALL Saxpby1LastNs(bnorm, xcstate, one, x0, pxc)
         CALL last_ntype_par
         CALL PadSides(pxc)
      END IF LACTIVE0

      CALL prepare_gmres_trial
      CALL funct3d_par(lscreen0, ier_flag_res)
      IF (lactive) CALL last_ns_par

      fsq_nl = fsqr + fsqz + fsql
      IF (l_coupled_gmres) fsq_nl = gmres_merit
      IF (trial_context) THEN
         IF (.NOT.valid_gmres_trial(.TRUE.))
     &      fsq_nl = SQRT(HUGE(one))
      END IF
      nfcn = nfcn + 1

      END SUBROUTINE GetNLForce_par
  !------------------------------------------------
  ! 
  !------------------------------------------------
      SUBROUTINE last_ns_par 
         USE xstuff

         REAL(dp), ALLOCATABLE, DIMENSION(:)  :: tmp
         ALLOCATE (tmp(ntmaxblocksize*ns))

         CALL tolastns(pgc,tmp)
         CALL copylastns(tmp,pgc)

         CALL tolastns(pxcdot,tmp)
         CALL copylastns(tmp,pxcdot)

         CALL tolastns(pxc,tmp)
         CALL copylastns(tmp,pxc)

         CALL tolastns(pxsave,tmp)
         CALL copylastns(tmp,pxsave)

         DEALLOCATE(tmp)

      END SUBROUTINE last_ns_par
  !------------------------------------------------

  !------------------------------------------------
      SUBROUTINE last_ntype_par 
         USE xstuff

         REAL(dp), ALLOCATABLE, DIMENSION(:)  :: tmp
         ALLOCATE (tmp(ntmaxblocksize*ns))

         CALL tolastntype(pgc,tmp)
         CALL copylastntype(tmp,pgc)

         CALL tolastntype(pxcdot,tmp)
         CALL copylastntype(tmp,pxcdot)

         CALL tolastntype(pxc,tmp)
         CALL copylastntype(tmp,pxc)

         CALL tolastntype(pxsave,tmp)
         CALL copylastntype(tmp,pxsave)

         DEALLOCATE(tmp)
      END SUBROUTINE last_ntype_par
  !------------------------------------------------

  !------------------------------------------------
      SUBROUTINE gmres_fun_par (ier_flag, itype, lscreen,
     &                          step_accepted, reuse_force)
      USE precon2d, ONLY: block_precond_par
      USE blocktridiagonalsolver, ONLY: l_colscale
      USE xstuff
      USE vmec_main, ONLY: ftolv
      USE vmec_params, ONLY: bad_jacobian_flag, input_error_flag
      USE gmres_lib, ONLY: gmres_par, gmres_info
      INTEGER, INTENT(IN) :: itype
      INTEGER, INTENT(OUT) :: ier_flag
      LOGICAL, INTENT(IN) :: lscreen
      LOGICAL, OPTIONAL, INTENT(OUT) :: step_accepted
      LOGICAL, OPTIONAL, INTENT(IN) :: reuse_force
      TYPE(gmres_info) :: gi
      INTEGER :: n, m, jlo, jhi, solver_status
      INTEGER :: icntl(9), info(3)
      REAL(dp) :: cntl(5), fact, fsq_base, fsq_trial
      REAL(dp) :: fsq_best
      LOGICAL :: ok, accepted, reuse_base

      ier_flag = 0
      IF (PRESENT(step_accepted)) step_accepted = .FALSE.
!     GMRESR and TFQMR operate only on the serial equilibrium arrays.
      IF (itype.NE.1) THEN
         ier_flag = input_error_flag
         RETURN
      END IF
      LGMRESCALL = .TRUE.
      accepted = .FALSE.
      n = neqs
      IF (lfirst) THEN
         lfirst = .FALSE.
         IF (grank.EQ.0) THEN
            WRITE(*,'(2x,a,/)') 'Beginning safeguarded GMRES iterations'
            WRITE(nthreed,'(2x,a,/)')
     &         'Beginning safeguarded GMRES iterations'
         END IF
      END IF
      CALL init_dgmres(icntl,cntl)
      cntl(1) = 1.e-3_dp
      icntl(2) = 0
      icntl(3) = 20
      icntl(4) = 0
      icntl(5) = 3
      icntl(6) = 0
      icntl(7) = 20
!     Sample inner iterates sparsely to seed the nonlinear line search.
      icntl(9) = 2
      info = 0
      gi%m=20; gi%icntl=icntl; gi%cntl=cntl; gi%info=info
      gi%ftol=ftolv
      ALLOCATE(gi%rcounts(nranks),gi%disp(nranks))
      gi%startglobrow=tlglob; gi%endglobrow=trglob
      gi%iam=rank; gi%nprocs=nranks
      gi%rcounts=ntblkrcounts; gi%disp=ntblkdisp
      gi%mblk_size=ntmaxblocksize
      gi%my_comm=NS_COMM; gi%my_comm_world=RUNVMEC_COMM_WORLD
      gi%lactive=lactive; gi%lverbose=lscreen

      CALL begin_gmres_trials
      CALL CopyLastNtype(pxc, pxsave)
      reuse_base = .FALSE.
      IF (PRESENT(reuse_force)) reuse_base = reuse_force
      CALL collective_trial_ok(reuse_base, .TRUE.)
      CALL prepare_gmres_trial
      IF (reuse_base) THEN
         fsq_base = gmres_merit
      ELSE
         CALL funct3d_par(lscreen0, ier_flag_res)
         nfcn = nfcn+1
         fsq_base = fsqr+fsqz+fsql
         IF (l_coupled_gmres) fsq_base = gmres_merit
      END IF
      IF (.NOT.valid_gmres_trial(.TRUE.)) GOTO 900
      IF (.NOT.ieee_is_finite(fsq_base)) GOTO 900
      CALL CopyLastNtype(pgc, pxcdot)
      CALL CopyLastNtype(pgc, pgc, -one)
      CALL last_ns_par
      IF (lactive) THEN
         jlo = (tlglob-1)*ntmaxblocksize+1
         jhi = trglob*ntmaxblocksize
      END IF
      CALL gmres_par(n, gi, matvec_par, block_precond_par,
     &               getnlforce_par, pxcdot, pgc)
      ok = .NOT.bad_matvec
      IF (lactive) ok = ok .AND.
     &                 ALL(ieee_is_finite(pxcdot(jlo:jhi)))
      CALL collective_trial_ok(ok, .TRUE.)
      solver_status = gi%info(1)
      CALL MPI_BCAST(solver_status, 1, MPI_INTEGER, 0,
     &               RUNVMEC_COMM_WORLD, MPI_ERR)
      CALL last_ntype_par
!     -4 is only the inner iteration limit: a finite inexact direction
!     may still be accepted by the physical nonlinear merit test.
      IF (.NOT.ok .OR. (solver_status.LT.0 .AND.
     &                 solver_status.NE.-4)) THEN
         IF (grank.EQ.0 .AND. lscreen)
     &      WRITE(*,*) 'Unusable GMRES direction, status:',
     &                 solver_status
         GOTO 900
      END IF

!     Accept the first valid, sufficiently decreasing physical residual.
!     This leaves all working arrays at the accepted point; do not scan
!     past it and then restore only scalar diagnostics from another point.
      fact = one
      fsq_best = HUGE(one)
      DO m=1,8
         CALL SaxpbyLastNtype(fact, pxcdot, one, pxsave, pxc)
         IF (lactive) CALL PadSides(pxc)
         CALL prepare_gmres_trial
         CALL funct3d_par(lscreen0, ier_flag_res)
         nfcn = nfcn+1
         ok = valid_gmres_trial(.TRUE.)
         IF (ok) THEN
            fsq_trial = fsqr+fsqz+fsql
            IF (l_coupled_gmres) fsq_trial = gmres_merit
            fsq_best = MIN(fsq_trial,fsq_best)
            ok = fsq_trial.LE.(one-1.e-4_dp*fact)*fsq_base
         END IF
         CALL collective_trial_ok(ok, .TRUE.)
         IF (ok) THEN
            accepted = .TRUE.
            EXIT
         END IF
         fact = 0.5_dp*fact
      END DO
      IF (grank.EQ.0 .AND. lscreen) THEN
         IF (accepted) THEN
            WRITE(*,930) fact, fsq_base, fsq_trial
         ELSE
            WRITE(*,931) fsq_base, fsq_best, solver_status
         END IF
      END IF
  930 FORMAT(2x,"GMRES step: alpha=",f8.5,2x,
     &       "merit before/after ",1p,2e13.5)
  931 FORMAT(2x,"GMRES rejected: base/best ",1p,2e13.5,
     &       2x,"linear status ",i4)
      IF (.NOT.accepted) GOTO 900
!     FUNCT3D evaluated x_base + D*(y-y_base). Store that physical x,
!     retaining its fresh force and diagnostics without a duplicate call.
      IF (ictrl_prec2d.EQ.1 .AND. l_colscale)
     &   CALL SaxLastNtype(pxcdot,pcol_scale,pxcdot)
      CALL SaxpbyLastNtype(fact, pxcdot, one, pxsave, pxc)
      CALL CopyLastNtype(pxc, pxsave)
      GOTO 950

 900  CONTINUE
!     Recompute the original point, since a rejected trial may have
!     overwritten every force/geometry/vacuum work array. The caller can
!     then resume the native relaxation from that point.
      CALL CopyLastNtype(pxsave, pxc)
      pxcdot = 0
      IF (lactive) CALL PadSides(pxc)
      CALL prepare_gmres_trial
      CALL funct3d_par(lscreen0, ier_flag_res)
      nfcn = nfcn+1
      ier_flag = ier_flag_res
      IF (.NOT.valid_gmres_trial(.TRUE.) .AND. ier_flag.EQ.0)
     &   ier_flag = bad_jacobian_flag
 950  CONTINUE
      CALL end_gmres_trials
      IF (PRESENT(step_accepted)) step_accepted = accepted
      DEALLOCATE(gi%rcounts,gi%disp)
      LGMRESCALL = .FALSE.
      END SUBROUTINE gmres_fun_par

      SUBROUTINE GetNLForce(xcstate, fsq_nl, bnorm)
      USE xstuff, ONLY: xc, gc, x0=>xsave
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp),INTENT(IN)  :: xcstate(neqs), bnorm
      REAL(dp),INTENT(OUT) :: fsq_nl
!-----------------------------------------------
!undo internal gmres normalization
      xc(1:neqs) = x0(1:neqs) + bnorm*xcstate(1:neqs)

      CALL prepare_gmres_trial
      CALL funct3d(lscreen0, ier_flag_res)
      fsq_nl = fsqr+fsqz+fsql
      IF (l_coupled_gmres) fsq_nl = gmres_merit
      IF (trial_context) THEN
         IF (.NOT.valid_gmres_trial(.FALSE.))
     &      fsq_nl = SQRT(HUGE(one))
      END IF

      nfcn = nfcn + 1

      END SUBROUTINE GetNLForce

      SUBROUTINE matvec (p, Ap, ndim)
      USE vmec_main, ONLY: irst
      USE vmec_params, ONLY: bad_jacobian_flag
      USE xstuff, ONLY: xc, x0=>xsave, gc0=>xcdot, gc
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)                    :: ndim
      REAL(dp), INTENT(in), DIMENSION(ndim)  :: p
      REAL(dp), INTENT(out), DIMENSION(ndim) :: Ap
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: l
      REAL(dp) :: delta, pmax
C-----------------------------------------------
!
!     Computes linearized matrix product A*p = [F(x0+delta*p) - F(x0)]/delta, about point x0
!     Must scale p so delta*max|p| ~ sqrt(epsilon) to get accurate numerical derivative
!
!      Ap = -p
!      GOTO 90

      delta = SQRT(EPSILON(delta))
      pmax = SQRT(SUM(p(1:ndim)**2))
      delta = delta/MAX(delta, pmax)

      xc(1:ndim) = x0(1:ndim) + delta*p(1:ndim)
      CALL prepare_gmres_trial
      CALL funct3d(lscreen0, ier_flag_res)
      Ap = (gc(1:ndim) - gc0(1:ndim))/delta

      IF (trial_context) THEN
         IF (.NOT.valid_gmres_trial(.FALSE.) .OR.
     &       .NOT.ALL(ieee_is_finite(Ap))) THEN
            bad_matvec = .TRUE.
            Ap = -p
         END IF
      ELSE IF (ier_flag_res.NE.0 .OR. irst.NE.1 .OR.
     &         .NOT.ALL(ieee_is_finite(Ap))) THEN
         IF (legacy_matvec_status.EQ.0) THEN
            legacy_matvec_status = ier_flag_res
            IF (legacy_matvec_status.EQ.0)
     &         legacy_matvec_status = bad_jacobian_flag
         END IF
         Ap = -p
      END IF

 90   CONTINUE
      nfcn = nfcn + 1

      END SUBROUTINE matvec

      SUBROUTINE gmres_fun (ier_flag, itype, step_accepted,
     &                      reuse_force)
      USE precon2d, ONLY: block_precond
      USE xstuff
      USE vmec_main, ONLY: ftolv
      USE vmec_params, ONLY: bad_jacobian_flag
      USE gmres_lib, ONLY: gmres_ser, gmres_info
      INTEGER, INTENT(IN) :: itype
      INTEGER, INTENT(OUT) :: ier_flag
      LOGICAL, OPTIONAL, INTENT(OUT) :: step_accepted
      LOGICAL, OPTIONAL, INTENT(IN) :: reuse_force
      TYPE(gmres_info) :: gi
      INTEGER :: n, m, icntl(9), info(3)
      REAL(dp) :: cntl(5), fact, fsq_base, fsq_trial, fsq_best
      LOGICAL :: ok, accepted, reuse_base

      ier_flag = 0
      IF (PRESENT(step_accepted)) step_accepted = .TRUE.
      IF (itype.EQ.2) THEN
         CALL gmresr_fun(ier_flag)
         RETURN
      ELSE IF (itype.EQ.3) THEN
         CALL qmr_fun(ier_flag)
         RETURN
      END IF
      accepted = .FALSE.
      IF (PRESENT(step_accepted)) step_accepted = .FALSE.
      n = neqs
      IF (lfirst) THEN
         lfirst = .FALSE.
         WRITE(*,'(2x,a,/)') 'Beginning safeguarded GMRES iterations'
         WRITE(nthreed,'(2x,a,/)')
     &      'Beginning safeguarded GMRES iterations'
      END IF
      CALL init_dgmres(icntl,cntl)
      cntl(1) = 1.e-3_dp
      icntl(2)=0; icntl(3)=20; icntl(4)=0; icntl(5)=3
      icntl(6)=0; icntl(7)=20
!     Sample inner iterates sparsely to seed the nonlinear line search.
      icntl(9) = 2
      info=0
      gi%m=20; gi%icntl=icntl; gi%cntl=cntl; gi%info=info
      gi%ftol=ftolv
      CALL begin_gmres_trials
      xsave = xc
      reuse_base = .FALSE.
      IF (PRESENT(reuse_force)) reuse_base = reuse_force
      CALL prepare_gmres_trial
      IF (reuse_base) THEN
         fsq_base = gmres_merit
      ELSE
         CALL funct3d(lscreen0, ier_flag_res)
         nfcn = nfcn+1
         fsq_base = fsqr+fsqz+fsql
         IF (l_coupled_gmres) fsq_base = gmres_merit
      END IF
      IF (.NOT.valid_gmres_trial(.FALSE.)) GOTO 900
      IF (.NOT.ieee_is_finite(fsq_base)) GOTO 900
      xcdot = gc
      gc = -gc
      IF (MAXVAL(ABS(gc)).EQ.0._dp) GOTO 900
      CALL gmres_ser(n, gi, matvec, block_precond, getnlforce,
     &               xcdot, gc)
      IF (bad_matvec .OR.
     &    (gi%info(1).LT.0 .AND. gi%info(1).NE.-4) .OR.
     &    .NOT.ALL(ieee_is_finite(xcdot))) GOTO 900
      fact = one
      fsq_best = HUGE(one)
      DO m=1,8
         xc(1:n) = xsave(1:n) + fact*xcdot(1:n)
         CALL prepare_gmres_trial
         CALL funct3d(lscreen0, ier_flag_res)
         nfcn = nfcn+1
         ok = valid_gmres_trial(.FALSE.)
         IF (ok) THEN
            fsq_trial = fsqr+fsqz+fsql
            IF (l_coupled_gmres) fsq_trial = gmres_merit
            fsq_best = MIN(fsq_trial,fsq_best)
            ok = fsq_trial.LE.(one-1.e-4_dp*fact)*fsq_base
         END IF
         IF (ok) THEN
            accepted = .TRUE.
            EXIT
         END IF
         fact = 0.5_dp*fact
      END DO
      IF (accepted) THEN
         WRITE(*,930) fact, fsq_base, fsq_trial
      ELSE
         WRITE(*,931) fsq_base, fsq_best, gi%info(1)
      END IF
  930 FORMAT(2x,"GMRES step: alpha=",f8.5,2x,
     &       "merit before/after ",1p,2e13.5)
  931 FORMAT(2x,"GMRES rejected: base/best ",1p,2e13.5,
     &       2x,"linear status ",i4)
      IF (.NOT.accepted) GOTO 900
      xsave = xc
      GOTO 950

 900  CONTINUE
      xc = xsave
      xcdot = 0
      CALL prepare_gmres_trial
      CALL funct3d(lscreen0, ier_flag_res)
      nfcn = nfcn+1
      ier_flag = ier_flag_res
      IF (.NOT.valid_gmres_trial(.FALSE.) .AND. ier_flag.EQ.0)
     &   ier_flag = bad_jacobian_flag
 950  CONTINUE
      CALL end_gmres_trials
      IF (PRESENT(step_accepted)) step_accepted = accepted
      END SUBROUTINE gmres_fun

      SUBROUTINE gmresr_fun (ier_flag)
      USE vmec_main, ONLY: lfreeb, irzloff
      USE vmec_params, ONLY: bad_jacobian_flag
      USE xstuff
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ier_flag
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ndim, jtrunc, mgmres, maxits, lwork, linear_status
      LOGICAL :: oktest
      REAL(dp) :: eps, resid
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: work, delx, brhs
      CHARACTER(len=3), PARAMETER :: stc="rel"
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     &                              'Beginning GMRESR iterations'
C-----------------------------------------------
      IF (lfirst) THEN
         lfirst = .false.
         WRITE (*,'(2x,a,/)') qmr_message
         WRITE (nthreed, '(2x,a,/)') qmr_message
      END IF

      oktest = .false.
      ndim   = neqs
      jtrunc = 10
      mgmres = 20
      maxits = 10
      lwork  = ndim*(2*jtrunc + mgmres + 2)
      eps = .3_dp

      ALLOCATE(work(lwork), delx(ndim), brhs(ndim), stat=ier_flag_res)
      IF (ier_flag_res .ne. 0) STOP 'Allocation failed in gmresr'

      xsave = xc
      xcdot = gc
      legacy_matvec_status = 0
      brhs = -gc(1:ndim)
      delx = 0

      CALL gmresr(oktest, ndim, jtrunc, mgmres, brhs, delx, work,
     &            eps, stc, maxits, resid, matvec, linear_status)

      ier_flag = legacy_matvec_status
      IF (.NOT.ALL(ieee_is_finite(delx)) .AND. ier_flag.EQ.0)
     &   ier_flag = bad_jacobian_flag
      xc = xsave
      IF (ier_flag.EQ.0) xc(1:ndim) = xsave(1:ndim) + delx(1:ndim)
      IF (.NOT.lfreeb) xc(ns:2*irzloff:ns)=xsave(ns:2*irzloff:ns)

      DEALLOCATE (work, delx, brhs)

      END SUBROUTINE gmresr_fun

      SUBROUTINE qmr_fun(ier_flag)
      USE vmec_dim, ONLY: ns, mpol1, ntor1
      USE vmec_params, ONLY: ntmax, bad_jacobian_flag
      USE vmec_main, ONLY: lfreeb, irzloff
      USE xstuff
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, INTENT(OUT) :: ier_flag
      INTEGER :: ndim, nlen, nlim, ierr, info(4)
      INTEGER :: revcom, colx, colb
      INTEGER :: nty, ntyp, mt, mp, nt, np, jp, js
      REAL(dp), DIMENSION(neqs,9) :: vecs
      REAL(dp) :: tol = 1.E-3_dp
      CHARACTER(LEN=*), PARAMETER :: qmr_message = 
     1                               'Beginning TF-QMR iterations'
      LOGICAL, PARAMETER :: ldump_fort33 = .false.
C-----------------------------------------------
      ndim = SIZE(vecs,1)
      nlen = ndim
      nlim = 10
      xcdot = gc
      xsave = xc
      legacy_matvec_status = 0

      IF (lfirst) THEN
         lfirst = .false.
         WRITE (*,'(2x,a,/)') qmr_message
         WRITE (nthreed, '(2x,a,/)') qmr_message
         IF (ldump_fort33) THEN

!     CHECK THAT INITIALLY, PRECONDITIONER YIELDS (APPROXIMATE) IDENTITY MATRIX
!     OUTER LOOP: SUM OF XC PERTURBATION
!     INNER LOOP: GETS ALLL THE FORCES IN RESPONSE TO EACH OUTER LOOP PERTURBATION
         ierr = 0
         DO nty = 1, 3*ntmax
            WRITE(33, '(a,i4)') "NTYPE' (XC-pert) = ",nty
         DO mt = 0, mpol1
            WRITE(33, '(a,i4)') "M' = ",mt
         DO nt = 0, ntor1-1
            WRITE(33, '(a,i4)') "N' = ",nt
         DO js = 1, ns
            ierr = ierr+1
            IF (ierr .gt. ndim) EXIT
            IF (MOD(ierr,50).eq.0) PRINT '(2x,a,f8.2,a)', 'Progress: ', 
     1                             REAL(100*ierr)/ndim, ' %'
            IF (js.eq.ns .and. .not.lfreeb) CYCLE
            colx = 1;  colb = 3
            vecs(:,colx) = 0; vecs(ierr,colx) = 1
            CALL matvec(vecs(1,colx), vecs(1,colb), ndim)
            WRITE (33, '(a,i4,2x,a,i5,2x,a,1p,e12.2)') "js' = ", js, 
     1        ' ipert = ',ierr,' Ap[ipert,ipert] = ', vecs(ierr, colb)
            colx = 0
            DO ntyp = 1, 3*ntmax
               DO mp = 0, mpol1
                  DO np = 0, ntor1-1
                     DO jp = 1, ns
                        colx = colx + 1
                        IF (colx .gt. ndim) CYCLE
                        IF (colx.eq.ierr .or. 
     1                     ABS(vecs(colx,colb)).lt.0.05_dp) CYCLE
                        WRITE (33, 123)'ntype = ', ntyp,' m = ',mp,
     1                  ' n = ', np,' js = ', jp,' iforce = ',colx,
     2                  ' Ap[iforce,ipert] = ', vecs(colx,colb)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         END DO
         END DO
         END DO
 
         PRINT '(/,2x,a,/)','Jacobian check in file FORT.33'

         END IF

      END IF
 123  FORMAT(4x,4(a,i4),a,i6,a,1p,e12.3)
!
!     INITIALIZE vecs
!
      vecs(:ndim,2) = -gc(:ndim)
      vecs(:ndim,3) =  gc(:ndim)

      ierr = 100000
      info = 0
      info(1) = ierr

 10   CALL dutfx (ndim,nlen,nlim,vecs,tol,info)
      revcom = info(2)
      colx   = info(3)
      colb   = info(4)
      IF (revcom .eq. 1) THEN
         CALL matvec (vecs(1,colx), vecs(1,colb), ndim)
         IF (legacy_matvec_status.EQ.0) GO TO 10
      END IF

      ier_flag = legacy_matvec_status
      IF (.NOT.ALL(ieee_is_finite(vecs(1:ndim,1))) .AND.
     &    ier_flag.EQ.0) ier_flag = bad_jacobian_flag
      xc = xsave
      IF (ier_flag.EQ.0)
     &   xc(1:ndim) = xsave(1:ndim) + vecs(1:ndim,1)
      IF (.NOT.lfreeb) xc(ns:2*irzloff:ns)=xsave(ns:2*irzloff:ns)

      END SUBROUTINE qmr_fun

      END MODULE gmres_mod

      SUBROUTINE Truncate(num, iprec0)
      USE stel_kinds, ONLY: dp
      IMPLICIT NONE
!     NEEDED TO RESOLVE CALL IN gmres_par
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)     :: iprec0
      REAL(dp), INTENT(INOUT) :: num
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER*24 :: chnum, tchnum
!-----------------------------------------------
!
!     TRUNCATES double-precision to precision iprec0 digits, keeping exponent range of double
!     WRITE TO INTERNAL FILE TO DO TRUNCATION
!
!UNCOMMENT THESE LINES TO ACTIVATE
!      WRITE (chnum, '(a,i2,a,i2,a)') '(1p,e',iprec0+7,'.',iprec0,')'
!      WRITE (tchnum, chnum) num

!      READ (tchnum, chnum) num

      END SUBROUTINE Truncate
