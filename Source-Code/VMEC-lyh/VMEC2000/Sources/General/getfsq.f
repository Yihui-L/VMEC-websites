      SUBROUTINE getfsq_par(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main, ONLY: ns, mnsize
      USE vmec_params, ONLY: ntmax
      USE parallel_include_module
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: medge
      REAL(dp), INTENT(OUT) :: gnormr, gnormz
      REAL(dp), INTENT(IN) :: gnorm
      REAL(dp), DIMENSION(mnsize,ns,ntmax), INTENT(IN) :: gcr, gcz
      REAL(dp) :: totals(2)

      IF (.NOT.lactive) RETURN
      CALL getfsq_local(gcr, gcz, totals, medge)
      CALL sum_surface_scalars(totals, SIZE(totals))
      gnormr = gnorm*totals(1)
      gnormz = gnorm*totals(2)
      END SUBROUTINE getfsq_par

      SUBROUTINE getfsq_local(gcr, gcz, totals, medge)
      USE vmec_main, ONLY: ns, ns1, mnsize
      USE vmec_params, ONLY: ntmax
      USE parallel_include_module, ONLY: dp, tlglob, trglob
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: medge
      REAL(dp), DIMENSION(mnsize,ns,ntmax), INTENT(IN) :: gcr, gcz
      REAL(dp), INTENT(OUT) :: totals(2)
      INTEGER :: l

      totals = 0
      DO l = tlglob, MIN(trglob, ns1+medge)
         totals(1) = totals(1) + SUM(gcr(:,l,:)**2)
         totals(2) = totals(2) + SUM(gcz(:,l,:)**2)
      END DO
      END SUBROUTINE getfsq_local

      SUBROUTINE sum_surface_scalars(values, nvalues)
      USE parallel_include_module
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nvalues
      REAL(dp), INTENT(INOUT) :: values(nvalues)
#if defined(MPI_OPT)
      REAL(dp) :: ton, toff

      IF (nranks.LE.1) RETURN
      CALL second0(ton)
      CALL MPI_Allreduce(MPI_IN_PLACE, values, nvalues, MPI_REAL8,
     &                   MPI_SUM, NS_COMM, MPI_ERR)
      CALL second0(toff)
      allreduce_time = allreduce_time + (toff-ton)
#endif
      END SUBROUTINE sum_surface_scalars

      SUBROUTINE getfsq(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main, ONLY: rprec, ns, ns1, mnsize
      USE vmec_params, ONLY: ntmax
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: medge
      REAL(dp), INTENT(out) :: gnormr, gnormz
      REAL(dp), INTENT(in)  :: gnorm
      REAL(dp), DIMENSION(ns,mnsize*ntmax), INTENT(in) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: jsmax
!-----------------------------------------------
      jsmax = ns1 + medge
      gnormr = gnorm * SUM(gcr(:jsmax,:)**2)
      gnormz = gnorm * SUM(gcz(:jsmax,:)**2)

      END SUBROUTINE getfsq
