      SUBROUTINE scalpot(bvec, ivacskip, info)
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
      USE vacmod
      USE hessian_lapack, ONLY: hessian_getrf
      USE parallel_include_module
      USE timer_sub
      USE vmec_input, ONLY: l_cpu_accel
!$    USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: ivacskip
      INTEGER, INTENT(out) :: info
      REAL(dp), INTENT(out) :: bvec(mnpd2)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: ip, istore, istart, istore_max, ndim
      INTEGER :: first_ip, last_ip, worker, workers
      REAL(dp), ALLOCATABLE :: grpmn(:), gstore(:), amatrix(:)
      REAL(dp), ALLOCATABLE :: green(:,:), greenp(:,:), work(:,:,:)
      REAL(dp) :: ton, toff, tonscal


C-----------------------------------------------
      CALL second0(tonscal)
      info = 0

      IF (.NOT.ALLOCATED(amatsav)) THEN
         STOP 'AMATSAV: Allocation error in scalpot'
      END IF

      ALLOCATE (grpmn(nuv3*mnpd2), stat=ip)
      IF (ip .NE. 0) STOP 'GRPMN: Allocation error in scalpot'

!
!     COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL
!     ON EXIT, BVEC CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE
!     AND GRPMN CONTAINS TRANSFORM OF NORMAL DERIVATIVE
!     OF THE GREENS FUNCTION [PKM, EQ.(2.15)]
!
!     FOR ivacskip != 0, USE PREVIOUSLY COMPUTED bvecsav FOR SPEED
!

      ndim = mnpd2/mnpd
      CALL analyt (grpmn, bvec, ivacskip, ndim)

      IF (ivacskip .NE. 0) THEN
         IF (vacuum_lu_valid) THEN
            bvec = bvec + bvecsav
         ELSE
            info = -1
         END IF
      ELSE

         istore_max = MIN(64,nuv3)
         workers = 1
!$       IF (l_cpu_accel) workers = MIN(omp_get_max_threads(),
!$   &      istore_max, MAX(1,nuv3max-nuv3min+1))

         ALLOCATE (green(nuv,MERGE(istore_max,1,workers.GT.1)),
     &             gstore(nuv),
     &             greenp(nuv,istore_max), work(nuv,2,workers),
     &             amatrix(mnpd2*mnpd2), stat=ip)
         IF (ip .NE. 0) STOP 'Allocation error in scalpot'
         gstore = 0
         istart = nuv3min - 1

!        The primed points in a tile write independent Green columns.
!        Keep the integral accumulation in increasing ip order.
         DO first_ip = nuv3min, nuv3max, istore_max
            last_ip = MIN(nuv3max,first_ip+istore_max-1)
            CALL second0(ton)
            IF (workers.EQ.1) THEN
               DO ip = first_ip, last_ip
                  istore = ip-first_ip+1
                  CALL greenf(green(:,1), greenp(:,istore), ip,
     &                        work(:,1,1), work(:,2,1))
                  gstore = gstore + bexni(ip)*green(:,1)
               END DO
            ELSE
!$OMP PARALLEL DO NUM_THREADS(MIN(workers,last_ip-first_ip+1))
!$OMP& PRIVATE(istore,worker) SCHEDULE(STATIC)
               DO ip = first_ip, last_ip
                  istore = ip-first_ip+1
                  worker = 1
!$                worker = omp_get_thread_num()+1
                  CALL greenf(green(:,istore), greenp(:,istore), ip,
     &                        work(:,1,worker), work(:,2,worker))
               END DO
!$OMP END PARALLEL DO

               DO ip = first_ip, last_ip
                  istore = ip-first_ip+1
                  gstore = gstore + bexni(ip)*green(:,istore)
               END DO
            END IF
            CALL second0(toff)
            timer_vac(tgreenf) = timer_vac(tgreenf) + (toff-ton)
            greenf_time = timer_vac(tgreenf)
            CALL fourp(grpmn, greenp, last_ip-first_ip+1,
     &                 istart, last_ip, ndim)
         END DO

         CALL second0(ton)
         IF (vlactive) THEN
            CALL MPI_Allreduce(MPI_IN_PLACE, gstore, SIZE(gstore),
     &                         MPI_REAL8, MPI_SUM, VAC_COMM, MPI_ERR)
         END IF
         CALL second0(toff)
         allreduce_time = allreduce_time + (toff - ton)
         timer_vac(tallr) = timer_vac(tallr) + (toff-ton)
!
!        COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
!        AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
!
         CALL fouri (grpmn, gstore, amatrix, amatsav, bvec,
     &               bvecsav, ndim)
         DEALLOCATE (green, greenp, gstore, work, amatrix)

         info = -1
         IF (ALL(ieee_is_finite(amatsav))) THEN
            CALL hessian_getrf(mnpd2, amatsav, ipiv, info)
         END IF
         vacuum_lu_valid = info.EQ.0

      END IF

      DEALLOCATE (grpmn, stat=ip)

!
!     FINAL REDUCTION OVER VAC_COMM DONE ONCE HERE
!
      CALL second0(ton)
      IF (vlactive) THEN
         CALL MPI_Allreduce(MPI_IN_PLACE, bvec, SIZE(bvec), MPI_REAL8,
     &                      MPI_SUM, VAC_COMM, MPI_ERR)
      END IF

      CALL second0(toff)
      allreduce_time = allreduce_time + (toff - ton)
      timer_vac(tanar) = timer_vac(tanar) + (toff-ton)

      scalpot_time = scalpot_time + (toff - tonscal)

      END SUBROUTINE scalpot
