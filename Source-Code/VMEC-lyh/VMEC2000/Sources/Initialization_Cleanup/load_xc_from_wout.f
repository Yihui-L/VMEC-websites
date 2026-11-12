      SUBROUTINE load_xc_from_wout(rmn, zmn, lmn, lreset,
     1    ntor_in, mpol1_in, ns, reset_file, ierr)
      USE read_wout_mod, ONLY: w_rmnc=>rmnc, w_zmns=>zmns, w_lmns=>lmns,
     1    w_rmns=>rmns, w_zmnc=>zmnc, w_lmnc=>lmnc, w_xm=>xm, w_xn=>xn,
     2    w_ns=>ns, w_ntor=>ntor, w_mpol=>mpol, w_nfp=>nfp,
     3    w_mnmax=>mnmax, w_lasym=>lasym, read_wout_file,
     4    read_wout_deallocate
      USE vmec_params, ONLY: mscale, nscale, ntmax, lamscale,
     1                       rcc, rss, rsc, rcs, zsc, zcs, zcc, zss
      USE vparams, ONLY: one, zero, rprec
      USE vmec_input, ONLY: lasym, nfp
      USE vmec_main, ONLY: lthreed, p5 => cp5, phipf
      USE realspace, ONLY: shalf, sqrts
      USE parallel_include_module, ONLY: grank, gnranks,
     1                                 RUNVMEC_COMM_WORLD, MPI_ERR
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
#if defined(MPI_OPT)
      USE mpi_inc
#endif
#if defined(NETCDF)
      USE netcdf, ONLY: nf90_open, nf90_close, nf90_nowrite,
     1  nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension,
     2  nf90_get_var, nf90_noerr, nf90_enotvar, nf90_einval
#endif
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ns, mpol1_in, ntor_in
      REAL(rprec), DIMENSION(ns,0:ntor_in,0:mpol1_in,ntmax),
     1   INTENT(out) :: rmn, zmn, lmn
      LOGICAL, INTENT(inout) :: lreset
      CHARACTER(LEN=*), INTENT(inout) :: reset_file
      INTEGER, INTENT(out) :: ierr
      INTEGER :: mn, m, n, n1, js, mnmax
      REAL(rprec) :: t1, t2, smj, spjm
      REAL(rprec), ALLOCATABLE :: seed(:,:,:), xm(:), xn(:), temp(:,:)
      LOGICAL :: lnc, exists
      CHARACTER(LEN=LEN_TRIM(reset_file)+10) :: filename
      LOGICAL :: used(-ntor_in:ntor_in,0:mpol1_in)
#if defined(NETCDF)
      INTEGER :: file_id
#endif

!     Read and reconstruct the entire seed on rank zero. Lambda's inverse
!     half-mesh interpolation is a recurrence from the magnetic axis.
      ierr = 0
      IF (grank.EQ.0) THEN
         INQUIRE(FILE=reset_file, EXIST=exists)
         IF (exists) THEN
            filename = reset_file
            n1 = LEN_TRIM(filename)
            lnc = .FALSE.
            IF (n1.GE.3) lnc = filename(n1-2:n1).EQ.'.nc'
         ELSE
            filename = 'wout'
            CALL parse_extension(filename, reset_file, lnc)
         END IF

         mnmax = 1 + ntor_in + mpol1_in*(2*ntor_in+1)
         ALLOCATE(xm(mnmax), xn(mnmax),
     1            seed(mnmax,ns,MERGE(6,3,lasym)))
         seed = zero
         IF (lnc) THEN
#if defined(NETCDF)
            CALL read_seed_netcdf
#else
            ierr = 1
#endif
         ELSE
            CALL read_seed_text
         END IF
         IF (ierr.EQ.0) CALL reconstruct_seed
         IF (ierr.NE.0) THEN
            WRITE(*,'(2A)') 'Unable to read compatible WOUT seed: ',
     1           TRIM(filename)
            WRITE(*,'(A,4(1X,I0))') 'Expected NS, MPOL, NTOR, NFP:',
     1           ns, mpol1_in+1, ntor_in, nfp
         END IF
      END IF

      reset_file = ' '
#if defined(MPI_OPT)
      IF (gnranks.GT.1) THEN
         CALL MPI_Bcast(ierr,1,MPI_INTEGER,0,RUNVMEC_COMM_WORLD,MPI_ERR)
         IF (ierr.EQ.0) THEN
            CALL MPI_Bcast(rmn,SIZE(rmn),MPI_REAL8,0,
     1                     RUNVMEC_COMM_WORLD,MPI_ERR)
            CALL MPI_Bcast(zmn,SIZE(zmn),MPI_REAL8,0,
     1                     RUNVMEC_COMM_WORLD,MPI_ERR)
            CALL MPI_Bcast(lmn,SIZE(lmn),MPI_REAL8,0,
     1                     RUNVMEC_COMM_WORLD,MPI_ERR)
         END IF
      END IF
#endif
      IF (ierr.EQ.0) lreset = .FALSE.

      CONTAINS

      SUBROUTINE reconstruct_seed
      ierr = 1
      IF (nfp.LE.0 .OR. .NOT.ieee_is_finite(lamscale)) RETURN
      IF (lamscale.LE.zero) RETURN
      IF (.NOT.ALL(ieee_is_finite(seed))) RETURN
      IF (.NOT.ALL(ieee_is_finite(xm))) RETURN
      IF (.NOT.ALL(ieee_is_finite(xn))) RETURN
      IF (ANY(xm.LT.zero) .OR. ANY(xm.GT.mpol1_in)) RETURN
      IF (ANY(ABS(xn).GT.REAL(ntor_in,rprec)*nfp)) RETURN

      used = .FALSE.
      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn))/nfp
         IF (xm(mn).NE.m .OR. xn(mn).NE.n*nfp) RETURN
         IF (m.EQ.0 .AND. n.LT.0) RETURN
         IF (used(n,m)) RETURN
         used(n,m) = .TRUE.
      END DO

      rmn = zero
      zmn = zero
      lmn = zero
      DO mn = 1, mnmax
         m = NINT(xm(mn))
         n = NINT(xn(mn))/nfp
         n1 = ABS(n)
         t1 = one/(mscale(m)*nscale(n1))
         t2 = t1
         IF (n .lt. 0) t2 = -t2
         IF (n .eq. 0) t2 = zero
         rmn(:ns,n1,m,rcc) = rmn(:ns,n1,m,rcc) + t1*seed(mn,:ns,1)
         zmn(:ns,n1,m,zsc) = zmn(:ns,n1,m,zsc) + t1*seed(mn,:ns,2)
         lmn(:ns,n1,m,zsc) = lmn(:ns,n1,m,zsc) + t1*seed(mn,:ns,3)
         IF (lthreed) THEN
         rmn(:ns,n1,m,rss) = rmn(:ns,n1,m,rss) + t2*seed(mn,:ns,1)
         zmn(:ns,n1,m,zcs) = zmn(:ns,n1,m,zcs) - t2*seed(mn,:ns,2)
         lmn(:ns,n1,m,zcs) = lmn(:ns,n1,m,zcs) - t2*seed(mn,:ns,3)
         END IF
         IF (lasym) THEN
         rmn(:ns,n1,m,rsc) = rmn(:ns,n1,m,rsc) + t1*seed(mn,:ns,4)
         zmn(:ns,n1,m,zcc) = zmn(:ns,n1,m,zcc) + t1*seed(mn,:ns,5)
         lmn(:ns,n1,m,zcc) = lmn(:ns,n1,m,zcc) + t1*seed(mn,:ns,6)
         IF (lthreed) THEN
         rmn(:ns,n1,m,rcs) = rmn(:ns,n1,m,rcs) - t2*seed(mn,:ns,4)
         zmn(:ns,n1,m,zss) = zmn(:ns,n1,m,zss) + t2*seed(mn,:ns,5)
         lmn(:ns,n1,m,zss) = lmn(:ns,n1,m,zss) + t2*seed(mn,:ns,6)
         END IF
         END IF
         IF (m .eq. 0) THEN
            zmn(:ns,n1,m,zsc) = zero
            lmn(:ns,n1,m,zsc) = zero
            IF (lthreed) rmn(:ns,n1,m,rss) = zero
            IF (lasym) THEN
            rmn(:ns,n1,m,rsc) = zero
            IF (lthreed) THEN
            zmn(:ns,n1,m,zss) = zero
            lmn(:ns,n1,m,zss) = zero
            END IF
            END IF
         END IF
      END DO


!     Undo the same m=1 representation change used in WOUT conversion.
      IF (mpol1_in.GE.1 .AND. (lthreed .OR. lasym))
     1   ALLOCATE(temp(ns,0:ntor_in))
      IF (mpol1_in.GE.1 .AND. lthreed) THEN
         temp = rmn(:,:,1,rss)
         rmn(:,:,1,rss) = p5*(temp + zmn(:,:,1,zcs))
         zmn(:,:,1,zcs) = p5*(temp - zmn(:,:,1,zcs))
      END IF
      IF (mpol1_in.GE.1 .AND. lasym) THEN
         temp = rmn(:,:,1,rsc)
         rmn(:,:,1,rsc) = p5*(temp + zmn(:,:,1,zcc))
         zmn(:,:,1,zcc) = p5*(temp - zmn(:,:,1,zcc))
      END IF

      lmn(1,:,0,:) = lmn(2,:,0,:)
      spjm = shalf(2)/sqrts(2)
      IF (mpol1_in.GE.1) lmn(1,:,1,:) = lmn(2,:,1,:)/spjm
      lmn(1,:,2:,:) = zero
      DO js = 2, ns
         smj = shalf(js)/sqrts(js)
         IF (js.GT.2) spjm = shalf(js)/sqrts(js-1)
         DO m = 0, mpol1_in, 2
            lmn(js,:,m,:) = 2*lmn(js,:,m,:) - lmn(js-1,:,m,:)
         END DO
         DO m = 1, mpol1_in, 2
            lmn(js,:,m,:) =
     1         (2*lmn(js,:,m,:) - spjm*lmn(js-1,:,m,:))/smj
         END DO
      END DO
!     WROUT multiplies internal lambda by lamscale/phipf before
!     interpolating onto the half mesh. Reverse that order here.
      DO js = 1, ns
         lmn(js,:,:,:) = (phipf(js)/lamscale)*lmn(js,:,:,:)
      END DO
      IF (ALL(ieee_is_finite(rmn)) .AND.
     1    ALL(ieee_is_finite(zmn)) .AND.
     2    ALL(ieee_is_finite(lmn))) ierr = 0
      END SUBROUTINE reconstruct_seed

      SUBROUTINE read_seed_text
      INTEGER :: unit
      OPEN(NEWUNIT=unit,FILE=filename,STATUS='old',ACTION='read',
     1     IOSTAT=ierr)
      IF (ierr.NE.0) RETURN
      CALL read_wout_file(unit, ierr)
      CLOSE(unit)
      IF (ierr.NE.0) RETURN
      IF (w_ns.NE.ns .OR. w_ntor.NE.ntor_in .OR.
     1    w_mpol.NE.mpol1_in+1 .OR. w_nfp.NE.nfp .OR.
     2    w_mnmax.NE.mnmax) THEN
         ierr = 1
      ELSE
         xm = w_xm
         xn = w_xn
         seed(:,:,1) = w_rmnc
         seed(:,:,2) = w_zmns
         seed(:,:,3) = w_lmns
         IF (lasym .AND. w_lasym) THEN
            seed(:,:,4) = w_rmns
            seed(:,:,5) = w_zmnc
            seed(:,:,6) = w_lmnc
         END IF
      END IF
      CALL read_wout_deallocate
      END SUBROUTINE read_seed_text

#if defined(NETCDF)
      SUBROUTINE read_seed_netcdf
      INTEGER :: k, id, values(6), close_status, asym_vars
      CHARACTER(LEN=8), PARAMETER :: scalar_names(6) =
     1   (/'ns      ','mpol    ','ntor    ','nfp     ',
     2     'mnmax   ','ier_flag'/)
      CHARACTER(LEN=4), PARAMETER :: spectra(6) =
     1   (/'rmnc','zmns','lmns','rmns','zmnc','lmnc'/)

!     Read only data needed for the initial guess, never the seed's
!     boundary flags, coils, or equilibrium profiles.
      ierr = nf90_open(filename,nf90_nowrite,file_id)
      IF (ierr.NE.nf90_noerr) RETURN
      DO k = 1, 6
         ierr = seed_variable(TRIM(scalar_names(k)),[INTEGER::],id)
         IF (ierr.NE.nf90_noerr) GOTO 900
         ierr = nf90_get_var(file_id,id,values(k))
         IF (ierr.NE.nf90_noerr) GOTO 900
      END DO
      IF (ANY(values(1:5).NE.[ns,mpol1_in+1,ntor_in,nfp,mnmax])
     1    .OR. .NOT.ANY(values(6).EQ.[0,2,4])) THEN
         ierr = nf90_einval
         GOTO 900
      END IF

      ierr = seed_variable('xm',[mnmax],id)
      IF (ierr.NE.nf90_noerr) GOTO 900
      ierr = nf90_get_var(file_id,id,xm)
      IF (ierr.NE.nf90_noerr) GOTO 900
      ierr = seed_variable('xn',[mnmax],id)
      IF (ierr.NE.nf90_noerr) GOTO 900
      ierr = nf90_get_var(file_id,id,xn)
      IF (ierr.NE.nf90_noerr) GOTO 900

      asym_vars = 0
      DO k = 1, SIZE(seed,3)
         ierr = seed_variable(spectra(k),[mnmax,ns],id)
         IF (k.GT.3 .AND. ierr.EQ.nf90_enotvar) CYCLE
         IF (ierr.NE.nf90_noerr) GOTO 900
         IF (k.GT.3) asym_vars = asym_vars + 1
         ierr = nf90_get_var(file_id,id,seed(:,:,k))
         IF (ierr.NE.nf90_noerr) GOTO 900
      END DO
      ierr = nf90_noerr
      IF (asym_vars.NE.0 .AND. asym_vars.NE.3) ierr = nf90_einval

 900  CONTINUE
      close_status = nf90_close(file_id)
      IF (ierr.EQ.nf90_noerr) ierr = close_status
      END SUBROUTINE read_seed_netcdf

      INTEGER FUNCTION seed_variable(name, shape, id) RESULT(status)
      CHARACTER(LEN=*), INTENT(in) :: name
      INTEGER, INTENT(in) :: shape(:)
      INTEGER, INTENT(out) :: id
      INTEGER :: ndims, dimids(2), dimlen, j

      status = nf90_inq_varid(file_id,name,id)
      IF (status.NE.nf90_noerr) RETURN
      status = nf90_inquire_variable(file_id,id,ndims=ndims)
      IF (status.NE.nf90_noerr) RETURN
      IF (ndims.NE.SIZE(shape)) THEN
         status = nf90_einval
         RETURN
      END IF
      IF (ndims.EQ.0) RETURN
      status = nf90_inquire_variable(file_id,id,dimids=dimids)
      IF (status.NE.nf90_noerr) RETURN
      DO j = 1, ndims
         status = nf90_inquire_dimension(file_id,dimids(j),len=dimlen)
         IF (status.NE.nf90_noerr) RETURN
         IF (dimlen.NE.shape(j)) THEN
            status = nf90_einval
            RETURN
         END IF
      END DO
      END FUNCTION seed_variable
#endif

      END SUBROUTINE load_xc_from_wout
