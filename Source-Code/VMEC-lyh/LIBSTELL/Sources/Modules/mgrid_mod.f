      MODULE mgrid_mod
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
      USE v3_utilities
      USE stel_kinds
#if defined(MPI_OPT)
      USE mpi_inc, ONLY: MPI_GROUP_NULL
#endif
      USE vmec_input, ONLY: nbfld, nflxs, lfreeb, lrecon
      USE vsvd0, ONLY: nigroup, nparts, npfcoil, nbcoilsp, nfloops,
     1                 nbctotp
      IMPLICIT NONE
      PRIVATE :: ieee_is_finite

      INTEGER, PARAMETER :: nlimset = 2       !number of different limiters
      CHARACTER(LEN=*), PARAMETER :: 
     1   vn_br0 = 'br', vn_bp0 = 'bp', vn_bz0 = 'bz',
     2   vn_ar0 = 'ar', vn_ap0 = 'ap', vn_az0 = 'az',
     3   vn_ir = 'ir', vn_jz = 'jz',
     4   vn_kp = 'kp', vn_nfp = 'nfp',
     5   vn_rmin='rmin', vn_rmax='rmax', vn_zmin='zmin', 
     6   vn_zmax='zmax', vn_coilgrp='coil_group' 
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_nextcur = 'nextcur',  vn_mgmode='mgrid_mode', 
     2  vn_coilcur = 'raw_coil_cur',
     3  vn_flp = 'nobser', vn_nobd = 'nobd', vn_nbset = 'nbsets', 
     4  vn_nbfld = 'nbfld',
     2  ln_flp = 'flux loops', ln_nobd = 'Connected flux loops',
     3  ln_nbset = 'B-coil loops', ln_next = 'External currents',
     4  ln_nbfld = 'B-coil measurements'
 
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!     nr0b, np0b, nz0b:
!              :  grid dimensions of magnetic field in mgrid file
!     nbvac    :  total number of grid points (nr0b*np0b*nz0b) in mgrid file
!     bvac(:,1):  br (radial component of external magnetic field)
!     bvac(:,2):  bp (toroidal component)
!     bvac(:,3) = bz (z-component)
!     rminb, rmaxb : min (max) radial dimension of grid in mgrid
!     zminb, zmaxb : min (max) vertical dimension of grid in mgrid
!
!     nextcur:         no. of EXTERNAL current groups (eg., TF, PF, helical)
!     raw_coil_current  array of raw currents for each coil group
!     mgrid_mode     = 'S', scaled mode; = 'R', raw mode
!     curlabel:   array of labels describing each current group
!                     included in green''s FUNCTION BFIELD response
!
!               - - - - - - - - - - - - - - - - - -
!               FOR DIAGNOSTICS AND DATA ANALYSIS
!               (HERE,COILS ARE FOR MEASURING FIELDS, FLUXES)
!    iconnect:   two-dimensional array describing electrical
!     needflx:    =NEEDIT, loop required for flux match
!                    >=ISYMCOIL, loop required but flux computed by
!                               invoking symmetry in Z
!                    =IDONTNEED, loop not required for flux match
!    needbfld:    =NEEDIT, loop required for B-field match
!                    =ISAMECOIL, loop at same position as previous loop
!                    >=ISYMCOIL, loop required but B-field computed
!                               by invoking symmetry in Z
!                    =IDONTNEED, loop not required in B-field match
!      dsiext:    connected flux loop signals due to EXTERNAL coils
!      plflux:    array of measured (inferred) plasma contrib. to flux loops
!      plbfld:    array of measured (inferred) plasma contrib. to B-loops
!                      connection of up to four flux loops. Specifies
!                     the sign and flux loop number of (up to) four
!                     connected individual loops (indexing based on
!                     xobser,zobser arrays).
!        nobd:   number of connected flux loop measurements
!      nobser:   number of individual flux loop positions
!      nbsets:   number of B-coil sets defined in mgrid file
!  nbcoils(n):   number of bfield coils in each set defined in mgrid file
!    nbcoilsn:   total number of bfield coils defined in mgrid file
!      nbfldn:   total number of EXTERNAL bfield measurements used in matching
!  bloopnames:   array of labels describing b-field sets
!    dsilabel:   array of labels describing connected flux loops
!      xobser:   array of flux loop R-positions
!      zobser:   array of flux loop Z-positions
! rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
! zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
! abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
!
      INTEGER :: nr0b, np0b, nfper0, nz0b
      INTEGER :: nobd, nobser, nextcur, nbfldn, nbsets, nbcoilsn
      INTEGER :: nbvac, nbcoil_max, nlim, nlim_max, nsets,
     1           nrgrid, nzgrid
      INTEGER, DIMENSION(:), POINTER :: needflx=>NULL(), nbcoils=>NULL()
      INTEGER, DIMENSION(:), POINTER :: limitr=>NULL(), nsetsn=>NULL()
      INTEGER, DIMENSION(:,:), POINTER ::
     1                     iconnect=>NULL(), needbfld=>NULL()
      REAL(rprec) :: rminb, zminb, rmaxb, zmaxb, delrb, delzb
      REAL(rprec) ::rx1, rx2, zy1, zy2, condif
      REAL(rprec), DIMENSION(:,:), POINTER :: bvac=>NULL()
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1                brvac=>NULL(), bzvac=>NULL(), bpvac=>NULL()
      LOGICAL, PRIVATE :: bvac_is_shared = .FALSE.
#if defined(MPI_OPT)
      REAL(rprec), DIMENSION(:,:), POINTER, PRIVATE ::
     1                        bvac_shared=>NULL()
      INTEGER, PRIVATE :: mgrid_group_old = MPI_GROUP_NULL
#endif
      REAL(rprec), DIMENSION(:,:), POINTER ::
     1   unpsiext=>NULL(), plbfld=>NULL(), rbcoil=>NULL(),
     2   zbcoil=>NULL(), abcoil=>NULL(),
     3   bcoil=>NULL(), rbcoilsqr=>NULL()
      REAL(rprec), DIMENSION(:), POINTER :: raw_coil_current=>NULL()
      REAL(rprec), DIMENSION(:), POINTER ::
     1   xobser=>NULL(), zobser=>NULL(), xobsqr=>NULL(),
     2   dsiext=>NULL(), psiext=>NULL(), plflux=>NULL(), b_chi=>NULL()
      CHARACTER(LEN=300) :: mgrid_path
      CHARACTER(LEN=300) :: mgrid_path_old = " "
!     BVAC is already weighted by EXTCUR and sampled on NV planes.
      LOGICAL, PRIVATE :: mgrid_cache_valid = .FALSE.
      INTEGER, PRIVATE :: mgrid_nv_old = -1, mgrid_nfp_old = -1
      INTEGER, PRIVATE :: mgrid_nextcur_old = 0
      INTEGER, PRIVATE :: mgrid_nbcoil_max_old = 0,
     1                    mgrid_nlim_max_old = 0
      INTEGER, PRIVATE :: mgrid_nflxs_old = 0
      INTEGER, PRIVATE :: mgrid_nbfld_old(SIZE(nbfld)) = 0
      LOGICAL, PRIVATE :: mgrid_lfreeb_old = .FALSE.
      CHARACTER(LEN=1), PRIVATE :: mgrid_mode_old = " "
      REAL(rprec), ALLOCATABLE, PRIVATE :: mgrid_extcur_old(:)
      INTEGER, PARAMETER, PRIVATE :: mgrid_ik = SELECTED_INT_KIND(18)
      CHARACTER(LEN=4096), PRIVATE :: mgrid_identity_old = " "
      INTEGER(mgrid_ik), PRIVATE :: mgrid_stamp_old(8) = 0
      PRIVATE :: mgrid_file_identity, mgrid_sync_error,
     1           mgrid_check_header, read_bin_diagnostics,
     2           allocate_bin_observations, free_mgrid_diagnostics,
     3           mgrid_invalidate_cache, mgrid_prepare_bvac
#if defined(MPI_OPT)
      PRIVATE :: mgrid_window_group, mgrid_release_shared_bvac,
     1           mgrid_cache_same_group
#endif
      CHARACTER(LEN=30), DIMENSION(:), POINTER :: curlabel=>NULL()
      CHARACTER(LEN=15), DIMENSION(:), POINTER :: 
     1                         dsilabel=>NULL(), bloopnames=>NULL()
      CHARACTER(LEN=30) :: tokid
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1                         dbcoil=>NULL(), pfcspec=>NULL()
      REAL(rprec), DIMENSION(:,:), POINTER :: 
     1    rlim=>NULL(), zlim=>NULL(), reslim=>NULL(), seplim=>NULL()
      CHARACTER(LEN=1) :: mgrid_mode

      INTEGER :: win_needflx, win_nbcoils, win_limitr, win_nsetsn,
     1           win_iconnect, win_needbfld, win_bvac,
     2           win_brvac, win_bzvac, win_bpvac

#if defined(NETCDF)
      PRIVATE :: read_mgrid_bin, read_mgrid_nc, mgrid_nc_values
#else
      PRIVATE :: read_mgrid_bin
#endif
      
      CONTAINS

      SUBROUTINE read_mgrid (mgrid_file, extcur, nv, nfp, lscreen, 
     1                       ier_flag, comm)
      USE system_mod
      USE mpi_inc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   V a r i a b l e s
C-----------------------------------------------
!     
!     mgrid_file:     full path to mgrid file
!     lscreen   :     logical controlling output to screen
!     ier_flag  ;     error flag returned to caller
!     extcur(n)    :  external current multiplier for bfield(n) components
!     comm      : Optional mpi communicator.
!
      INTEGER, INTENT(out)          :: ier_flag
      INTEGER, INTENT(in)           :: nv, nfp
      LOGICAL, INTENT(in)           :: lscreen
      REAL(rprec), INTENT(in)       :: extcur(:)
      CHARACTER(len=*), INTENT(in)  :: mgrid_file
      INTEGER, INTENT(in), OPTIONAL :: comm
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
#if defined(VMS)
      CHARACTER(LEN=*), PARAMETER :: mgrid_defarea='vmec$:[makegrid]'
#else
      CHARACTER(LEN=*), PARAMETER :: mgrid_defarea='$HOME/vmec/MAKEGRID'
#endif
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
!
!     lgrid_exist  :   logical set if mgrid file is found in given path
!
      INTEGER :: istat, ii
      CHARACTER(LEN=200) :: home_dir
      LOGICAL :: lgrid_exist, lfind, lmpi_local, lmpi_init
      LOGICAL :: cache_hit, all_cache_hit
      INTEGER :: comm_local
      CHARACTER(LEN=4096) :: file_identity, after_identity
      INTEGER(mgrid_ik) :: file_stamp(8), after_stamp(8)
      LOGICAL :: identity_known, after_known
C-----------------------------------------------

      ier_flag = 0
      lmpi_local = .FALSE.
      comm_local = 0
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         comm_local = comm
         CALL MPI_INITIALIZED(lmpi_init, istat)
         lmpi_local = lmpi_init
      END IF
#endif

      IF (nv.LE.0 .OR. nfp.LE.0) ier_flag = 9
      IF (.NOT.ALL(ieee_is_finite(extcur))) ier_flag = 9
      CALL mgrid_sync_error(ier_flag, comm)
      IF (ier_flag.NE.0) RETURN

      mgrid_path = TRIM(mgrid_file)

      INQUIRE (file=mgrid_path,exist=lgrid_exist,iostat=istat)
      IF (istat.ne.0 .or. .not.lgrid_exist) THEN
          IF (lscreen) PRINT *,' MGRID FILE NOT FOUND IN SPECIFIED ',
     1       'PATH: SEARCHING DEFAULT AREA'
          ii = INDEX(mgrid_file,'/',back=.true.)
          istat = INDEX(mgrid_defarea, '$HOME')
          IF (istat .ne. 0) THEN
             CALL getenv('HOME', home_dir)
             IF (istat .gt. 1) THEN
                home_dir = mgrid_defarea(1:istat-1) // TRIM(home_dir)
     1                   // mgrid_defarea(istat+5:)
             ELSE
                home_dir = TRIM(home_dir) // mgrid_defarea(istat+5:)
             END IF
          ELSE
             home_dir = mgrid_defarea
          END IF
#if defined(VMS)
          mgrid_path = TRIM(home_dir) // mgrid_file(ii+1:)
#else
          mgrid_path = TRIM(home_dir) // '/' // mgrid_file(ii+1:)
#endif
          INQUIRE (file=mgrid_path,exist=lgrid_exist,iostat=istat)
      END IF

!     Missing/unreadable paths must not make only some ranks skip the
!     collective reader or silently change the requested boundary model.
      IF (istat.NE.0 .OR. .NOT.lgrid_exist) ier_flag = 9
      CALL mgrid_sync_error(ier_flag, comm)
      IF (ier_flag.NE.0) THEN
         CALL mgrid_invalidate_cache
         RETURN
      END IF
      CALL mgrid_file_identity(mgrid_path, file_identity, file_stamp,
     1                         identity_known)

!     A cache hit is valid only for the same complete, successful read.
!     Keep the allocation/shape tests separate: Fortran need not
!     short-circuit logical expressions.
      cache_hit = .FALSE.
      IF (mgrid_cache_valid .AND. identity_known .AND.
     1    file_identity.EQ.mgrid_identity_old .AND.
     2    mgrid_path.EQ.mgrid_path_old .AND.
     3    ALL(file_stamp.EQ.mgrid_stamp_old) .AND.
     3    nv.EQ.mgrid_nv_old .AND. nfp.EQ.mgrid_nfp_old .AND.
     4    ASSOCIATED(curlabel) .AND. ASSOCIATED(bvac) .AND.
     5    nflxs.EQ.mgrid_nflxs_old .AND.
     6    ALL(nbfld.EQ.mgrid_nbfld_old) .AND.
     7    (lfreeb.EQV.mgrid_lfreeb_old) .AND.
     8    (lmpi_local.EQV.bvac_is_shared)) THEN
         IF (ALLOCATED(mgrid_extcur_old)) THEN
            IF (SIZE(extcur).EQ.SIZE(mgrid_extcur_old)) THEN
               cache_hit = ALL(extcur.EQ.mgrid_extcur_old)
            END IF
         END IF
      END IF
#if defined(MPI_OPT)
!     Every rank must take the same branch before collective MGRID I/O.
      IF (lmpi_local) THEN
         IF (cache_hit) cache_hit = mgrid_cache_same_group(comm_local)
         CALL MPI_ALLREDUCE(cache_hit, all_cache_hit, 1, MPI_LOGICAL,
     1                      MPI_LAND, comm_local, istat)
         IF (istat.NE.MPI_SUCCESS) THEN
            ier_flag = 9
            RETURN
         END IF
         cache_hit = all_cache_hit
      END IF
#endif
      IF (cache_hit) THEN
!        VSETUP resets these public values before a new VMEC run.
         nextcur = mgrid_nextcur_old
         mgrid_mode = mgrid_mode_old
         nbcoil_max = mgrid_nbcoil_max_old
         nlim_max = mgrid_nlim_max_old
         IF (lscreen) PRINT *,' mgrid file previously parsed!'
         RETURN
      END IF
!     Do not leave a valid key pointing to partially overwritten data.
      CALL mgrid_invalidate_cache
      CALL free_mgrid_diagnostics

      IF (lgrid_exist) THEN
         IF (lscreen) PRINT '(2x,2a)',
     1     'Opening vacuum field file: ', TRIM(mgrid_file)
!
!        Parse mgrid file name, look for .nc extension (netcdf format)
! 
         ii = LEN_TRIM(mgrid_path) - 2
         lfind = .FALSE.
         IF (ii.GE.1) lfind = (mgrid_path(ii:ii+2) == '.nc')
         IF (lfind) THEN
#if defined(NETCDF)
            IF (lmpi_local) THEN
               CALL read_mgrid_nc (mgrid_path, extcur, nv, nfp, 
     1                          ier_flag, lscreen, comm_local)
            ELSE
               CALL read_mgrid_nc (mgrid_path, extcur, nv, nfp, 
     1                          ier_flag, lscreen)
            END IF
#else
            ier_flag = 9
#endif
         ELSE
            IF (lmpi_local) THEN
               CALL read_mgrid_bin (mgrid_path, extcur, nv, nfp,
     1                          ier_flag, lscreen, comm_local)
            ELSE
               CALL read_mgrid_bin (mgrid_path, extcur, nv, nfp,
     1                          ier_flag, lscreen)
            END IF
         END IF

         CALL mgrid_sync_error(ier_flag, comm)
         IF (ier_flag.NE.0) RETURN
!SPH060517         IF (np0b .ne. nv) THEN
         IF (MOD(np0b, nv).NE.0) THEN
            PRINT *,' NZETA=',nv,
     1      ' DOES NOT DIVIDE EVENLY INTO NP0B=',np0b,' IN MGRID FILE'
            ier_flag = 9
         ELSE IF (nfper0.ne.nfp) THEN
            PRINT *,' NFP(READ in) = ',nfp,' DOES NOT AGREE WITH ',
     1      'NFPER (in vacuum field file) = ',nfper0
            ier_flag = 9
         END IF

      END IF
      
      IF (ier_flag .ne. 0) RETURN

      IF (.not.lgrid_exist .or. ier_flag.ne.0) THEN
         lfreeb = .false.
         lrecon = .false.
         IF (lscreen) THEN
            PRINT *, ' Error opening/reading mgrid file in dir: ',
     1                TRIM(home_dir)
            PRINT *, ' User must supply vacuum bfield in mgrid to ',
     1                'run vmec in free-boundary mode!'
            PRINT *, ' Proceeding to run vmec in',
     1               ' fixed boundary mode'
         END IF
      ELSE
         IF (ALLOCATED(mgrid_extcur_old))
     1      DEALLOCATE(mgrid_extcur_old)
         ALLOCATE(mgrid_extcur_old(SIZE(extcur)))
         mgrid_extcur_old = extcur
         mgrid_nv_old = nv
         mgrid_nfp_old = nfp
         mgrid_nextcur_old = nextcur
         mgrid_nbcoil_max_old = nbcoil_max
         mgrid_nlim_max_old = nlim_max
         mgrid_nflxs_old = nflxs
         mgrid_nbfld_old = nbfld
         mgrid_lfreeb_old = lfreeb
         mgrid_mode_old = mgrid_mode
         mgrid_path_old = mgrid_path
!        Do not cache if the platform cannot identify the file.
         CALL mgrid_file_identity(mgrid_path, after_identity,
     1                            after_stamp, after_known)
         mgrid_identity_old = file_identity
         mgrid_stamp_old = file_stamp
         mgrid_cache_valid = identity_known .AND. after_known .AND.
     1      after_identity.EQ.file_identity .AND.
     2      ALL(after_stamp.EQ.file_stamp)
#if defined(MPI_OPT)
         IF (mgrid_cache_valid .AND. lmpi_local) THEN
            CALL MPI_COMM_GROUP(comm_local, mgrid_group_old, istat)
            IF (istat.NE.MPI_SUCCESS) CALL mgrid_invalidate_cache
         END IF
#endif
      END IF

      END SUBROUTINE read_mgrid

      SUBROUTINE mgrid_invalidate_cache
#if defined(MPI_OPT)
      USE mpi_inc
      INTEGER :: ierr
      IF (mgrid_group_old.NE.MPI_GROUP_NULL) THEN
         CALL MPI_GROUP_FREE(mgrid_group_old,ierr)
      END IF
#endif
      mgrid_cache_valid = .FALSE.
      mgrid_path_old = ' '
      END SUBROUTINE mgrid_invalidate_cache

#if defined(MPI_OPT)
      LOGICAL FUNCTION mgrid_cache_same_group(comm)
      USE mpi_inc
      INTEGER, INTENT(in) :: comm
      INTEGER :: caller_group, relation, ierr
      mgrid_cache_same_group = .FALSE.
      IF (mgrid_group_old.EQ.MPI_GROUP_NULL) RETURN
      CALL MPI_COMM_GROUP(comm,caller_group,ierr)
      IF (ierr.NE.MPI_SUCCESS) RETURN
      CALL MPI_GROUP_COMPARE(caller_group,mgrid_group_old,relation,ierr)
      IF (ierr.EQ.MPI_SUCCESS) mgrid_cache_same_group =
     1                          relation.EQ.MPI_IDENT
      CALL MPI_GROUP_FREE(caller_group,ierr)
      END FUNCTION mgrid_cache_same_group

!     A shared window may be released only with its whole owning group.
!     Check locally before the caller combines errors across its ranks.
      SUBROUTINE mgrid_window_group(comm, identical, complete, ier)
      USE mpi_inc
      INTEGER, INTENT(in) :: comm
      LOGICAL, INTENT(out) :: identical, complete
      INTEGER, INTENT(inout) :: ier
      INTEGER :: caller_group, window_group, common_group, relation
      INTEGER :: ierr, common_size, window_size

      identical = .FALSE.
      complete = .TRUE.
      IF (.NOT.ASSOCIATED(bvac_shared)) RETURN
      CALL MPI_COMM_GROUP(comm,caller_group,ierr)
      IF (ierr.NE.MPI_SUCCESS) THEN
         ier = 9
         RETURN
      END IF
      CALL MPI_WIN_GET_GROUP(win_bvac,window_group,ierr)
      IF (ierr.NE.MPI_SUCCESS) THEN
         CALL MPI_GROUP_FREE(caller_group,ierr)
         ier = 9
         RETURN
      END IF
      CALL MPI_GROUP_COMPARE(caller_group,window_group,relation,ierr)
      IF (ierr.NE.MPI_SUCCESS) THEN
         ier = 9
      ELSE IF (relation.EQ.MPI_IDENT) THEN
         identical = .TRUE.
      ELSE IF (relation.EQ.MPI_UNEQUAL) THEN
         CALL MPI_GROUP_INTERSECTION(caller_group,window_group,
     1                               common_group,ierr)
         IF (ierr.EQ.MPI_SUCCESS) THEN
            CALL MPI_GROUP_SIZE(common_group,common_size,ierr)
            CALL MPI_GROUP_SIZE(window_group,window_size,ierr)
            complete = common_size.EQ.window_size
            CALL MPI_GROUP_FREE(common_group,ierr)
         ELSE
            ier = 9
         END IF
      END IF
      CALL MPI_GROUP_FREE(caller_group,ierr)
      CALL MPI_GROUP_FREE(window_group,ierr)
      END SUBROUTINE mgrid_window_group

      SUBROUTINE mgrid_release_shared_bvac(ier)
      USE mpi_inc
      INTEGER, INTENT(inout) :: ier
      INTEGER :: ierr

      IF (.NOT.ASSOCIATED(bvac_shared)) RETURN
      CALL MPI_WIN_FREE(win_bvac,ierr)
      IF (ierr.NE.MPI_SUCCESS) THEN
         ier = 9
         RETURN
      END IF
      NULLIFY(bvac_shared)
      IF (bvac_is_shared) THEN
         NULLIFY(bvac,brvac,bpvac,bzvac)
         bvac_is_shared = .FALSE.
      END IF
      END SUBROUTINE mgrid_release_shared_bvac
#endif

!     Keep a window separate from the active view. A serial caller can
!     switch to private storage without a collective on the old group.
      SUBROUTINE mgrid_prepare_bvac(npoints, shar_rank, shar_comm,
     1                              shared, ier, comm)
      USE mpi_sharmem
      USE mpi_inc
      INTEGER, INTENT(in) :: npoints, shar_rank, shar_comm
      LOGICAL, INTENT(in) :: shared
      INTEGER, INTENT(inout) :: ier
      INTEGER, INTENT(in), OPTIONAL :: comm
      INTEGER :: istat
#if defined(MPI_OPT)
      LOGICAL :: same_group, complete, reusable, all_reusable
#endif

      NULLIFY(brvac,bpvac,bzvac)
      istat = 0
#if defined(MPI_OPT)
      IF (shared) THEN
         CALL mgrid_window_group(shar_comm,same_group,complete,ier)
         IF (.NOT.complete) ier = 9
         CALL mgrid_sync_error(ier,comm)
         IF (ier.NE.0) RETURN

         reusable = .FALSE.
         IF (same_group) reusable = SIZE(bvac_shared,1).EQ.npoints
         CALL MPI_ALLREDUCE(reusable,all_reusable,1,MPI_LOGICAL,
     1                      MPI_LAND,shar_comm,istat)
         IF (istat.NE.MPI_SUCCESS) ier = 9
         CALL mgrid_sync_error(ier,comm)
         IF (ier.NE.0) RETURN
         IF (.NOT.all_reusable) CALL mgrid_release_shared_bvac(ier)

         IF (ASSOCIATED(bvac) .AND. .NOT.bvac_is_shared) THEN
            DEALLOCATE(bvac,STAT=istat)
            IF (istat.NE.0) ier = 9
         END IF
         CALL mgrid_sync_error(ier,comm)
         IF (ier.NE.0) RETURN

         IF (.NOT.ASSOCIATED(bvac_shared)) THEN
            CALL MPIALLOC(bvac_shared,npoints,3,shar_rank,0,
     1                    shar_comm,win_bvac)
         END IF
         IF (.NOT.ASSOCIATED(bvac_shared)) ier = 9
         CALL mgrid_sync_error(ier,comm)
         IF (ier.NE.0) RETURN
         bvac => bvac_shared
         bvac_is_shared = .TRUE.
         RETURN
      END IF
#endif
      IF (bvac_is_shared) NULLIFY(bvac)
      bvac_is_shared = .FALSE.
      IF (ASSOCIATED(bvac)) THEN
         IF (SIZE(bvac,1).NE.npoints) DEALLOCATE(bvac,STAT=istat)
      END IF
      IF (istat.EQ.0 .AND. .NOT.ASSOCIATED(bvac))
     1   ALLOCATE(bvac(npoints,3),STAT=istat)
      IF (istat.NE.0) ier = 9
      CALL mgrid_sync_error(ier,comm)
      END SUBROUTINE mgrid_prepare_bvac

!     Resolve cwd aliases; identify replacement and modification before hits.
!     The supplied Linux profile uses the fixed Linux statx UAPI, including
!     64-bit inode/size and nanosecond mtime/ctime. Unknown identity means
!     a conservative cache miss, never a pathname-only cache hit.
      SUBROUTINE mgrid_file_identity(path, identity, stamp, known)
#if defined(VMEC_LINUX_STATX)
      USE, INTRINSIC :: ISO_C_BINDING
#endif
      CHARACTER(LEN=*), INTENT(in) :: path
      CHARACTER(LEN=*), INTENT(out) :: identity
      INTEGER(mgrid_ik), INTENT(out) :: stamp(8)
      LOGICAL, INTENT(out) :: known
      CHARACTER(LEN=4096) :: cwd
      INTEGER :: ierr
#if defined(VMEC_LINUX_STATX)
!     linux/stat.h: all offsets are fixed by the 256-byte kernel ABI.
      TYPE, BIND(C) :: statx_time
         INTEGER(C_INT64_T) :: seconds
         INTEGER(C_INT32_T) :: nanoseconds, reserved
      END TYPE
      TYPE, BIND(C) :: statx_data
         INTEGER(C_INT32_T) :: mask, blksize
         INTEGER(C_INT64_T) :: attributes
         INTEGER(C_INT32_T) :: nlink, uid, gid
         INTEGER(C_INT16_T) :: mode, spare0
         INTEGER(C_INT64_T) :: inode, bytes, blocks, attributes_mask
         TYPE(statx_time) :: atime, btime, ctime, mtime
         INTEGER(C_INT32_T) :: rdev_major, rdev_minor, dev_major,
     1                        dev_minor
         INTEGER(C_INT64_T) :: spare(14)
      END TYPE
      INTERFACE
         FUNCTION file_statx(dirfd, pathname, flags, mask, buffer)
     1      BIND(C, NAME='statx') RESULT(status)
            IMPORT :: C_INT, C_CHAR, statx_data
            INTEGER(C_INT), VALUE :: dirfd, flags, mask
            CHARACTER(C_CHAR), INTENT(in) :: pathname(*)
            TYPE(statx_data), INTENT(out) :: buffer
            INTEGER(C_INT) :: status
         END FUNCTION
      END INTERFACE
      TYPE(statx_data) :: info
      INTEGER(C_INT) :: c_status
      INTEGER(C_INT), PARAMETER :: required = 960 ! 0x3c0
#endif
      identity = ' '
      stamp = 0
      known = .FALSE.
#if defined(__GFORTRAN__) || defined(VMEC_GNU_FILE_IDENTITY)
      IF (LEN_TRIM(path).EQ.0) RETURN
      IF (path(1:1).EQ.'/') THEN
         IF (LEN_TRIM(path).GT.LEN(identity)) RETURN
         identity = TRIM(path)
      ELSE
         CALL GETCWD(cwd, ierr)
         IF (ierr.NE.0) RETURN
         IF (LEN_TRIM(cwd)+1+LEN_TRIM(path).GT.LEN(identity)) RETURN
         identity = TRIM(cwd)//'/'//TRIM(path)
      END IF
#if defined(VMEC_LINUX_STATX)
      IF (C_SIZEOF(info).NE.256) RETURN
      c_status=file_statx(-100_C_INT,TRIM(identity)//C_NULL_CHAR,
     1                    0_C_INT,required,info)
      IF (c_status.NE.0) RETURN
      IF (IAND(info%mask,required).NE.required) RETURN
      stamp=(/INT(info%dev_major,mgrid_ik),
     1        INT(info%dev_minor,mgrid_ik),info%inode,info%bytes,
     2        info%mtime%seconds,INT(info%mtime%nanoseconds,mgrid_ik),
     3        info%ctime%seconds,INT(info%ctime%nanoseconds,mgrid_ik)/)
#else
!     GNU 11 STAT silently truncates large inode numbers. Platforms
!     without the opt-in 64-bit interface therefore reread conservatively.
      RETURN
#endif
      known = .TRUE.
#endif
      END SUBROUTINE mgrid_file_identity

!     All members must call this at the same reader stage.
      SUBROUTINE mgrid_sync_error(ier, comm)
      USE mpi_inc
      INTEGER, INTENT(inout) :: ier
      INTEGER, INTENT(in), OPTIONAL :: comm
      INTEGER :: combined, mgrid_mpi_status
      LOGICAL :: initialized
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_INITIALIZED(initialized, mgrid_mpi_status)
         IF (initialized) THEN
            CALL MPI_ALLREDUCE(ier, combined, 1, MPI_INTEGER,
     1                         MPI_MAX, comm, mgrid_mpi_status)
            ier = combined
            IF (mgrid_mpi_status.NE.MPI_SUCCESS) ier = 9
         END IF
      END IF
#endif
      END SUBROUTINE mgrid_sync_error

!     Reject invalid dimensions before divisions, indexing or collectives.
      SUBROUTINE mgrid_check_header(nv,nfp,ncur,ier,comm)
      USE mpi_inc
      INTEGER, INTENT(in) :: nv,nfp,ncur
      INTEGER, INTENT(inout) :: ier
      INTEGER, INTENT(in), OPTIONAL :: comm
      INTEGER :: dims(7),lo(7),hi(7),mgrid_mpi_status
      LOGICAL :: initialized
      IF (nr0b.LT.2 .OR. nz0b.LT.2 .OR. np0b.LT.1) ier=9
      IF (nextcur.LT.1 .OR. nextcur.GT.nigroup) ier=9
      IF (nextcur.GT.ncur .OR. nfper0.NE.nfp) ier=9
      IF (nv.LE.0) THEN
         ier=9
      ELSE
         IF (MOD(np0b,nv).NE.0) ier=9
      END IF
      IF (.NOT.ALL(ieee_is_finite((/rminb,rmaxb,zminb,zmaxb/)))) THEN
         ier=9
      ELSE IF (rmaxb.LE.rminb .OR. zmaxb.LE.zminb) THEN
         ier=9
      END IF
      IF (INT(nr0b,mgrid_ik)*INT(nz0b,mgrid_ik)*nv.GT.
     1    INT(HUGE(nbvac),mgrid_ik)) ier=9
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_INITIALIZED(initialized,mgrid_mpi_status)
         IF (initialized) THEN
            dims=(/nr0b,nz0b,np0b,nfper0,nextcur,nv,nfp/)
            CALL MPI_ALLREDUCE(dims,lo,7,MPI_INTEGER,MPI_MIN,
     1                         comm,mgrid_mpi_status)
            CALL MPI_ALLREDUCE(dims,hi,7,MPI_INTEGER,MPI_MAX,
     1                         comm,mgrid_mpi_status)
            IF (ANY(lo.NE.hi)) ier=9
         END IF
      END IF
#endif
      CALL mgrid_sync_error(ier,comm)
      END SUBROUTINE mgrid_check_header

      SUBROUTINE read_mgrid_python(mgrid_file, extcur, nv, nfp, nextcur)
      IMPLICIT NONE
      INTEGER, INTENT(in)           :: nv, nfp,nextcur
      REAL(rprec), INTENT(in)       :: extcur(1:nextcur)
      CHARACTER(len=*), INTENT(in)  :: mgrid_file
      INTEGER                       :: ier_flag
      ier_flag = 0
      CALL read_mgrid(mgrid_file,extcur,nv,nfp,.FALSE.,ier_flag)
      RETURN
      END SUBROUTINE read_mgrid_python

      
      SUBROUTINE read_mgrid_bin (filename, extcur, nv, nfp, ier_flag, 
     1                           lscreen,comm)
      USE, INTRINSIC :: iso_fortran_env, ONLY: iostat_end
      USE safe_open_mod
      USE mpi_sharmem
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y  A r g u m e n t s
C-----------------------------------------------
!
!     lstyle_2000  :   logical controlling ordering of magnetic field components read in
!
      INTEGER, INTENT(in)  :: nv, nfp
      CHARACTER(LEN=*), INTENT(in) :: filename
      REAL(rprec), INTENT(in) :: extcur(:)
      INTEGER, INTENT(inout)  :: ier_flag
      LOGICAL, INTENT(in)  :: lscreen
      INTEGER, INTENT(in), OPTIONAL :: comm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: 
     1           brtemp, bztemp, bptemp
      INTEGER :: iunit = 50
      INTEGER :: istat, ig, i, j, n, n1, nsets_max, k
      LOGICAL :: lstyle_2000, file_open
      LOGICAL :: lMPIInit, lraw_local, lmpi_local
      INTEGER :: mpi_rank, mpi_size, MPI_ERR
      INTEGER :: shar_rank, shar_comm
      INTEGER :: nskip, sh(1), diagnostic_dims(9)
      LOGICAL :: has_obs, has_limiter, end_of_tail
C-----------------------------------------------
      lraw_local = .false.; istat = 0
      file_open = .FALSE.
      lMPIInit = .FALSE.
      ier_flag = 0
      mpi_rank = 0; mpi_size = 1
      shar_comm = -1; shar_rank = 0
      lmpi_local = .FALSE.
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         CALL MPI_INITIALIZED(lMPIInit, MPI_ERR)
         IF (lMPIInit) THEN
            lmpi_local = .TRUE.
            CALL MPI_COMM_RANK(comm, mpi_rank, istat)
            CALL MPI_COMM_SIZE(comm, mpi_size, istat)
            CALL MPI_COMM_SPLIT_TYPE(comm, 
     1              MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, 
     2              shar_comm, istat)
            CALL MPI_COMM_RANK(shar_comm, shar_rank, istat)
         END IF
      END IF
#endif
!
!     OPEN FILE
!
      IF (shar_rank .eq. 0) THEN
         CALL safe_open(iunit, istat, filename,'old','unformatted')
         file_open = (istat.EQ.0)
      END IF
#if defined(MPI_OPT)
      IF (lmpi_local) 
     1   CALL MPI_BCAST(istat,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
#endif
      IF (istat.NE.0) ier_flag=9
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
!
!     READ IN BASIC DATA
!
      nr0b=0; nz0b=0; np0b=0; nfper0=0; nextcur=0
      rminb=0; rmaxb=0; zminb=0; zmaxb=0
      IF (shar_rank .eq. 0) THEN
         READ (iunit,iostat=istat) nr0b, nz0b, np0b, nfper0, nextcur
         IF (istat .ne. 0) ier_flag = 9
         READ(iunit,iostat=istat) rminb, zminb, rmaxb, zmaxb
         IF (istat .ne. 0) ier_flag = 9
      END IF
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
          CALL MPI_BCAST(nr0b,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(nz0b,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(np0b,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(nfper0,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(nextcur,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(rminb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(zminb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(rmaxb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(zmaxb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
          CALL MPI_BCAST(ier_flag,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
      END IF
#endif
!
!     CALCULATE HELPERS AND CHECK NFP AND NV
!
      lstyle_2000 = (nextcur.LT.0)
      IF (nextcur.EQ.-HUGE(nextcur)-1) THEN
         ier_flag=9
         nextcur=0
      ELSE
         nextcur=ABS(nextcur)
      END IF
      CALL mgrid_check_header(nv,nfp,SIZE(extcur),ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      delrb=(rmaxb-rminb)/(nr0b-1)
      delzb=(zmaxb-zminb)/(nz0b-1)
      nbvac=nr0b*nz0b*nv
!
!     ALLOCATE NEXTCUR SIZED ARRAYS
!
      istat = 0
      IF (.NOT. ASSOCIATED(curlabel)) THEN
         ALLOCATE (curlabel(nextcur), stat=istat)
      ELSE IF (SIZE(curlabel) .ne. nextcur) THEN
         DEALLOCATE (curlabel)
         ALLOCATE (curlabel(nextcur), stat=istat)
      END IF
!
!     READ NEXTCUR SIZED ARRAYS
!
      IF (istat.NE.0) ier_flag=9
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      curlabel=' '
      IF (shar_rank .eq. 0) THEN
         READ(iunit,iostat=istat) (curlabel(n),n=1,nextcur)
      END IF
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         CALL MPI_BCAST(istat,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(curlabel,SIZE(curlabel)*LEN(curlabel),
     1                  MPI_CHARACTER,
     1                  0,shar_comm,MPI_ERR)
!      IF (lraw_local)
!     1  CALL MPI_BCAST(raw_coil_current,nextcur,MPI_REAL8,0,
!     1                 shar_comm,MPI_ERR)
      END IF
#endif
!
!     ALLOCATE BVAC ARRAY
!
      IF (istat.NE.0) ier_flag=9
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      CALL mgrid_prepare_bvac(nbvac,shar_rank,shar_comm,
     1                        lmpi_local,ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
!
!     INITIALIZE BVAC ARRAY
!
      IF (shar_rank .eq. 0) bvac = 0
!
!     DEALLOCATE HELPER ARRAYS
!
      IF (ALLOCATED(brtemp)) DEALLOCATE(brtemp)
      IF (ALLOCATED(bptemp)) DEALLOCATE(bptemp)
      IF (ALLOCATED(bztemp)) DEALLOCATE(bztemp)
!
!     ALLOCATE BVAC HELPER ARRAYS
!
      nskip = np0b/nv
      IF (shar_rank .eq. 0) THEN
         ALLOCATE (brtemp(nr0b,nz0b,np0b), bptemp(nr0b,nz0b,np0b),
     1             bztemp(nr0b,nz0b,np0b), stat=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) THEN
         brtemp = 0
         bptemp = 0
         bztemp = 0
         END IF
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
!
!     READ BVAC HELPER ARRAYS
!
!     Only the shared-memory owner may initialize BVAC.  Synchronize
!     before summation so another rank cannot erase accumulated data.
      IF (shar_rank .eq. 0) bvac = 0
#if defined(MPI_OPT)
      IF (lmpi_local) CALL MPI_BARRIER(shar_comm,MPI_ERR)
#endif
      nskip = np0b/nv
      sh(1) = nbvac
      
      IF (shar_rank == 0) THEN
         DO ig = 1,nextcur
            IF (lstyle_2000) THEN
               READ(iunit, iostat=istat) brtemp, bptemp, bztemp
            ELSE
               READ(iunit, iostat=istat) (((brtemp(i,j,k),
     1                                   bztemp(i,j,k),bptemp(i,j,k), 
     2                               i= 1,nr0b),j=1,nz0b), k=1,np0b)
            END IF
            IF (istat.NE.0) THEN
               ier_flag=9
               EXIT
            END IF
            IF (.NOT.ALL(ieee_is_finite(brtemp)) .OR.
     1          .NOT.ALL(ieee_is_finite(bptemp)) .OR.
     2          .NOT.ALL(ieee_is_finite(bztemp))) THEN
               ier_flag=9
               EXIT
            END IF
!
!        STORE SUMMED BFIELD (OVER COIL GROUPS) IN BVAC
!         
         !CALL sum_bfield(bvac(1,1), brtemp, extcur(ig), nv)
         !CALL sum_bfield(bvac(1,2), bptemp, extcur(ig), nv)
         !CALL sum_bfield(bvac(1,3), bztemp, extcur(ig), nv)
         
            bvac(:,1) = bvac(:,1) + extcur(ig) * 
     1               RESHAPE(brtemp(:,:,1:np0b:nskip),sh)
         
            bvac(:,2) = bvac(:,2) + extcur(ig) * 
     1               RESHAPE(bptemp(:,:,1:np0b:nskip),sh)
         
            bvac(:,3) = bvac(:,3) + extcur(ig) * 
     1               RESHAPE(bztemp(:,:,1:np0b:nskip),sh)

         END DO
         IF (.NOT.ALL(ieee_is_finite(bvac))) ier_flag=9
      END IF
      np0b = nv

      IF (PRESENT(comm)) THEN
#if defined(MPI_OPT)
         CALL MPI_BCAST(istat,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
#endif
      END IF

      IF (istat.NE.0) ier_flag=9
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000

      IF (ALLOCATED(brtemp)) DEALLOCATE(brtemp)
      IF (ALLOCATED(bptemp)) DEALLOCATE(bptemp)
      IF (ALLOCATED(bztemp)) DEALLOCATE(bztemp)

      CALL assign_bptrs(bvac)

!     The mode/current records precede optional diagnostics in the
!     newer native layout. Only a clean EOF can omit the mode.
      mgrid_mode = 'N'
      end_of_tail = .FALSE.
      IF (lstyle_2000 .AND. shar_rank.EQ.0) THEN
         READ(iunit,IOSTAT=istat) mgrid_mode
         lraw_local = istat.EQ.0
         end_of_tail = istat.EQ.iostat_end
         IF (istat.NE.0 .AND. .NOT.end_of_tail) ier_flag = 9
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         CALL MPI_BCAST(mgrid_mode,1,MPI_CHARACTER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(lraw_local,1,MPI_LOGICAL,0,shar_comm,MPI_ERR)
      END IF
#endif
      IF (ASSOCIATED(raw_coil_current)) DEALLOCATE(raw_coil_current)
      IF (lraw_local) THEN
         ALLOCATE(raw_coil_current(nextcur),STAT=istat)
         IF (istat.NE.0) ier_flag = 9
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      IF (lraw_local .AND. shar_rank.EQ.0) THEN
         raw_coil_current = 0
         READ(iunit,IOSTAT=istat) raw_coil_current
         IF (istat.EQ.iostat_end) THEN
!           Mode-only legacy files fall back to the unknown mode.
            mgrid_mode = 'N'
            end_of_tail = .TRUE.
         ELSE IF (istat.NE.0) THEN
            ier_flag = 9
         END IF
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
#if defined(MPI_OPT)
      IF (lmpi_local .AND. lraw_local) THEN
         CALL MPI_BCAST(mgrid_mode,1,MPI_CHARACTER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(raw_coil_current,nextcur,MPI_REAL8,0,
     1                  shar_comm,MPI_ERR)
      END IF
#endif

!     Each node owner validates the full optional tail before peers
!     allocate buffers or enter payload collectives.
      has_obs = .FALSE.
      has_limiter = .FALSE.
      nsets_max = 0
      IF (shar_rank.EQ.0 .AND. .NOT.end_of_tail)
     1   CALL read_bin_diagnostics(iunit,has_obs,has_limiter,
     2                             nsets_max,ier_flag)
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         CALL MPI_BCAST(has_obs,1,MPI_LOGICAL,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(has_limiter,1,MPI_LOGICAL,0,shar_comm,MPI_ERR)
         diagnostic_dims = (/nobser,nobd,nbsets,nbcoil_max,nbfldn,
     1                       nlim,nlim_max,nsets,nsets_max/)
         CALL MPI_BCAST(diagnostic_dims,SIZE(diagnostic_dims),
     1                  MPI_INTEGER,0,shar_comm,MPI_ERR)
         IF (shar_rank.NE.0) THEN
            nobser=diagnostic_dims(1); nobd=diagnostic_dims(2)
            nbsets=diagnostic_dims(3); nbcoil_max=diagnostic_dims(4)
            nbfldn=diagnostic_dims(5); nlim=diagnostic_dims(6)
            nlim_max=diagnostic_dims(7); nsets=diagnostic_dims(8)
            nsets_max=diagnostic_dims(9)
            IF (has_obs) CALL allocate_bin_observations(ier_flag)
            IF (has_limiter .AND. ier_flag.EQ.0) THEN
               ALLOCATE(rlim(nlim_max,nlim),zlim(nlim_max,nlim),
     1                  reslim(nlim_max,nlim),seplim(nlim_max,nlim),
     2                  pfcspec(nparts,nsets_max,nsets),STAT=istat)
               IF (istat.NE.0) ier_flag=9
            END IF
         END IF
      END IF
#endif
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000

#if defined(MPI_OPT)
      IF (lmpi_local .AND. has_obs) THEN
         CALL MPI_BCAST(nbcoils,nbsets,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(xobser,nobser,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(zobser,nobser,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(dsilabel,SIZE(dsilabel)*LEN(dsilabel),
     1                  MPI_CHARACTER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(iconnect,SIZE(iconnect),MPI_INTEGER,0,
     1                  shar_comm,MPI_ERR)
         CALL MPI_BCAST(bloopnames,SIZE(bloopnames)*LEN(bloopnames),
     1                  MPI_CHARACTER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(rbcoil,SIZE(rbcoil),MPI_REAL8,0,
     1                  shar_comm,MPI_ERR)
         CALL MPI_BCAST(zbcoil,SIZE(zbcoil),MPI_REAL8,0,
     1                  shar_comm,MPI_ERR)
         CALL MPI_BCAST(abcoil,SIZE(abcoil),MPI_REAL8,0,
     1                  shar_comm,MPI_ERR)
         CALL MPI_BCAST(unpsiext,SIZE(unpsiext),MPI_REAL8,0,
     1                  shar_comm,MPI_ERR)
         CALL MPI_BCAST(dbcoil,SIZE(dbcoil),MPI_REAL8,0,
     1                  shar_comm,MPI_ERR)
      END IF
      IF (lmpi_local .AND. has_limiter) THEN
         CALL MPI_BCAST(rlim,SIZE(rlim),MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(zlim,SIZE(zlim),MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(pfcspec,SIZE(pfcspec),MPI_REAL8,0,
     1                  shar_comm,MPI_ERR)
         CALL MPI_BCAST(rx1,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(rx2,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(zy1,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(zy2,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(condif,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(nrgrid,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(nzgrid,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(tokid,30,MPI_CHARACTER,0,shar_comm,MPI_ERR)
      END IF
#endif
      IF (has_obs) THEN
         xobsqr = SQRT(xobser)
         nbcoilsn = SUM(nbcoils)
         DO n=1,nbsets
            rbcoilsqr(:nbcoils(n),n) = SQRT(rbcoil(:nbcoils(n),n))
         END DO
      ELSE IF (lscreen) THEN
         PRINT *, ' No observation data in mgrid file'
      END IF
      IF (has_limiter) THEN
         reslim = 0
         seplim = 0
      ELSE IF (has_obs .AND. lscreen) THEN
         PRINT *, ' No limiter data in mgrid file'
      END IF

      IF (file_open) CLOSE(iunit,IOSTAT=istat)
      file_open=.FALSE.

      

#if defined(MPI_OPT)
      IF (lmpi_local) THEN
        CALL MPI_COMM_FREE(shar_comm,istat)
      END IF
#endif

!
!     SUM UP CONTRIBUTIONS FROM INDIVIDUAL COIL GROUPS
!
      IF (lfreeb) THEN
         IF (nobser .gt. 0) psiext(:nobser) = 0
         IF (nbcoil_max.gt.0 .and. nbsets.gt.0)
     1       bcoil(:nbcoil_max, :nbsets) = 0

         DO ig = 1,nextcur
            IF (nobser .gt. 0)
     1      psiext(:nobser) = psiext(:nobser) +
     2                        extcur(ig)*unpsiext(:nobser,ig)
            DO n=1,nbsets
               n1 = nbcoils(n)
               bcoil(:n1,n) = bcoil(:n1,n) +
     1                        extcur(ig)*dbcoil(:n1,n,ig)
            ENDDO
         ENDDO
      ENDIF                   !!IF LFREEB

 9000 CONTINUE
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) THEN
         CALL free_mgrid_diagnostics
         IF (ASSOCIATED(raw_coil_current))
     1      DEALLOCATE(raw_coil_current)
         mgrid_mode = 'N'
      END IF
      IF (file_open) CLOSE(iunit,IOSTAT=istat)
      IF (ALLOCATED(brtemp)) DEALLOCATE(brtemp)
      IF (ALLOCATED(bptemp)) DEALLOCATE(bptemp)
      IF (ALLOCATED(bztemp)) DEALLOCATE(bztemp)
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         IF (shar_comm.NE.MPI_COMM_NULL)
     1      CALL MPI_COMM_FREE(shar_comm,istat)
      END IF
#endif
      END SUBROUTINE read_mgrid_bin

      SUBROUTINE allocate_bin_observations(ier_flag)
      INTEGER, INTENT(INOUT) :: ier_flag
      INTEGER :: istat
      ALLOCATE(nbcoils(nbsets), xobser(nobser), zobser(nobser),
     1       dsilabel(nobd), iconnect(4,nobser+nobd),
     2       unpsiext(nobser,nextcur), xobsqr(nobser), needflx(nobser),
     3       plflux(nobser+nobd), dsiext(nobd), psiext(nobser),
     4       bloopnames(nbsets), needbfld(nbcoil_max,nbsets),
     5       plbfld(nbcoil_max,nbsets), rbcoil(nbcoil_max,nbsets),
     6       zbcoil(nbcoil_max,nbsets), abcoil(nbcoil_max,nbsets),
     7       bcoil(nbcoil_max,nbsets), rbcoilsqr(nbcoil_max,nbsets),
     8       b_chi(nbsets), dbcoil(nbcoil_max,nbsets,nextcur),
     9       STAT=istat)
      IF (istat.NE.0) ier_flag = 9
      END SUBROUTINE allocate_bin_observations

      SUBROUTINE read_bin_diagnostics(iunit,has_obs,has_limiter,
     1                                nsets_max,ier_flag)
      USE, INTRINSIC :: iso_fortran_env, ONLY: iostat_end
      INTEGER, INTENT(IN) :: iunit
      INTEGER, INTENT(OUT) :: nsets_max, ier_flag
      LOGICAL, INTENT(OUT) :: has_obs, has_limiter
      INTEGER :: istat, n, n1, m, ig, i, j, k, dims(3)
      INTEGER :: coil_counts(SIZE(nbfld)), limiter_counts(nlimset)
      INTEGER :: pf_counts(nigroup)
      INTEGER(mgrid_ik) :: filesize, coord_count
      CHARACTER(LEN=15) :: label

      has_obs = .FALSE.
      has_limiter = .FALSE.
      nsets_max = 0
      ier_flag = 9
      dims = -1
      READ(iunit,IOSTAT=istat) dims
      IF (istat.EQ.iostat_end .AND. ALL(dims.EQ.-1)) THEN
         ier_flag = 0
         RETURN
      END IF
      IF (istat.NE.0) RETURN
      nobser = dims(1); nobd = dims(2); nbsets = dims(3)
      IF (MINVAL(dims).LT.0 .OR. nobser.GT.nfloops .OR.
     1    nobd.GT.nfloops .OR. nbsets.GT.SIZE(nbfld) .OR.
     2    nflxs.GT.nfloops) RETURN
      IF (ANY(nbfld(:nbsets).LT.0)) RETURN
      IF (SUM(INT(nbfld(:nbsets),mgrid_ik)).GT.nbctotp) RETURN
      nbfldn = SUM(nbfld(:nbsets))
      READ(iunit,IOSTAT=istat) coil_counts(:nbsets)
      IF (istat.NE.0) RETURN
      IF (ANY(coil_counts(:nbsets).LT.0) .OR.
     1    ANY(coil_counts(:nbsets).GT.nbcoilsp)) RETURN
      nbcoil_max = MAX(0,MAXVAL(coil_counts(:nbsets)))
      ier_flag = 0
      CALL allocate_bin_observations(ier_flag)
      IF (ier_flag.NE.0) RETURN
      ier_flag = 9
      nbcoils = coil_counts(:nbsets)
      iconnect = 0
      bloopnames = ' '
      rbcoil = 0; zbcoil = 0; abcoil = 0; dbcoil = 0

      READ(iunit,IOSTAT=istat) (xobser(n),zobser(n),n=1,nobser)
      IF (istat.NE.0) RETURN
      IF (ANY(.NOT.ieee_is_finite(xobser)) .OR.
     1    ANY(xobser.LT.0)) RETURN
      READ(iunit,IOSTAT=istat) dsilabel
      IF (istat.NE.0) RETURN
      READ(iunit,IOSTAT=istat) iconnect(:,:nobd)
      IF (istat.NE.0) RETURN
      DO n=1,nbsets
         IF (nbcoils(n).EQ.0) CYCLE
         READ(iunit,IOSTAT=istat) n1,label
         IF (istat.NE.0) RETURN
         IF (n1.LT.1 .OR. n1.GT.nbsets) RETURN
         bloopnames(n1) = label
         READ(iunit,IOSTAT=istat) (rbcoil(m,n),zbcoil(m,n),
     1                             abcoil(m,n),m=1,nbcoils(n))
         IF (istat.NE.0) RETURN
         IF (ANY(.NOT.ieee_is_finite(rbcoil(:nbcoils(n),n))) .OR.
     1       ANY(rbcoil(:nbcoils(n),n).LT.0)) RETURN
      END DO
      DO ig=1,nextcur
         READ(iunit,IOSTAT=istat) unpsiext(:,ig)
         IF (istat.NE.0) RETURN
         DO n=1,nbsets
            READ(iunit,IOSTAT=istat) dbcoil(:nbcoils(n),n,ig)
            IF (istat.NE.0) RETURN
         END DO
      END DO
      has_obs = .TRUE.

!     Read the bounded prefix first, then consume the whole record.
      nlim = -1
      READ(iunit,IOSTAT=istat) nlim
      IF (istat.EQ.iostat_end .AND. nlim.EQ.-1) THEN
         nlim = 0
         ier_flag = 0
         RETURN
      END IF
      IF (istat.NE.0) RETURN
      IF (nlim.LT.0 .OR. nlim.GT.nlimset) RETURN
      BACKSPACE(iunit,IOSTAT=istat)
      IF (istat.NE.0) RETURN
      READ(iunit,IOSTAT=istat) nlim,limiter_counts(:nlim)
      IF (istat.NE.0) RETURN
      IF (ANY(limiter_counts(:nlim).LT.0)) RETURN
      nlim_max = MAX(0,MAXVAL(limiter_counts(:nlim)))
      coord_count = 2*SUM(INT(limiter_counts(:nlim),mgrid_ik))
      INQUIRE(UNIT=iunit,SIZE=filesize,IOSTAT=istat)
      IF (istat.NE.0) RETURN
      IF (filesize.GE.0 .AND.
     1    coord_count.GT.filesize/(STORAGE_SIZE(rminb)/8)) RETURN
      IF (INT(nlim_max,mgrid_ik)*nlim.GT.HUGE(nlim_max)) RETURN
      ALLOCATE(rlim(nlim_max,nlim),zlim(nlim_max,nlim),
     1         reslim(nlim_max,nlim),seplim(nlim_max,nlim),STAT=istat)
      IF (istat.NE.0) RETURN
      rlim = 0; zlim = 0
      READ(iunit,IOSTAT=istat)
     1   ((rlim(i,j),zlim(i,j),i=1,limiter_counts(j)),j=1,nlim)
      IF (istat.NE.0) RETURN

      READ(iunit,IOSTAT=istat) nsets
      IF (istat.NE.0) RETURN
      IF (nsets.LT.0 .OR. nsets.GT.nigroup) RETURN
      BACKSPACE(iunit,IOSTAT=istat)
      IF (istat.NE.0) RETURN
      READ(iunit,IOSTAT=istat) nsets,pf_counts(:nsets)
      IF (istat.NE.0) RETURN
      IF (ANY(pf_counts(:nsets).LT.0) .OR.
     1    ANY(pf_counts(:nsets).GT.npfcoil)) RETURN
      nsets_max = MAX(0,MAXVAL(pf_counts(:nsets)))
      ALLOCATE(pfcspec(nparts,nsets_max,nsets),STAT=istat)
      IF (istat.NE.0) RETURN
      pfcspec = 0
      READ(iunit,IOSTAT=istat) (((pfcspec(i,j,k),i=1,nparts),
     1                              j=1,pf_counts(k)),k=1,nsets)
      IF (istat.NE.0) RETURN
      READ(iunit,IOSTAT=istat) rx1,rx2,zy1,zy2,condif,
     1                        nrgrid,nzgrid,tokid
      IF (istat.NE.0) RETURN
      has_limiter = .TRUE.
      ier_flag = 0
      END SUBROUTINE read_bin_diagnostics



#if defined(NETCDF)
!
!     PARALLEL MPI MODIFICATIONS ADDED BY Mark R. Cianciosa <cianciosamr@ornl.gov>, 092315
!
      SUBROUTINE read_mgrid_nc (filename, extcur, nv, nfp,
     1                          ier_flag, lscreen, comm)
      USE ezcdf
      USE mpi_inc
      USE mpi_sharmem
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(in) :: nv, nfp
      REAL(rprec), INTENT(in) :: extcur(:)
      INTEGER, INTENT(in), OPTIONAL :: comm
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(in) :: lscreen
      REAL(rprec), ALLOCATABLE :: brtemp(:,:,:), bptemp(:,:,:),
     1                           bztemp(:,:,:), bttemp(:,:)
      INTEGER :: ngrid, istat, ig, mpi_rank, mpi_size, MPI_ERR
      INTEGER :: shar_rank, shar_comm, temp_comm, temp_rank, temp_size
      INTEGER :: dimlens(3)
      CHARACTER(LEN=100) :: temp
      LOGICAL :: lmpi_local, lMPIInit, lraw_local, file_open

      ier_flag = 0
      istat = 0
      lmpi_local = .FALSE.
      lraw_local = .FALSE.
      file_open = .FALSE.
      mpi_rank = 0
      mpi_size = 1
      shar_rank = 0
      shar_comm = -1
      temp_comm = -1
      temp_rank = 0
      temp_size = 1
#if defined(MPI_OPT)
      shar_comm = MPI_COMM_NULL
      temp_comm = MPI_COMM_NULL
      IF (PRESENT(comm)) THEN
         CALL MPI_INITIALIZED(lMPIInit, MPI_ERR)
         IF (lMPIInit) THEN
            lmpi_local = .TRUE.
            CALL MPI_COMM_RANK(comm, mpi_rank, istat)
            CALL MPI_COMM_SIZE(comm, mpi_size, istat)
            CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0,
     1                              MPI_INFO_NULL, shar_comm, istat)
            CALL MPI_COMM_RANK(shar_comm, shar_rank, istat)
         END IF
      END IF
#endif
!     Each shared-memory root opens the metadata file. No rank leaves
!     a stage until errors have been combined over the full caller group.
      IF (shar_rank.EQ.0) THEN
         CALL cdf_open(ngrid, filename, 'r', istat)
         file_open = (istat.EQ.0)
         IF (istat.NE.0) ier_flag = 9
      END IF
      CALL mgrid_sync_error(ier_flag, comm)
      IF (ier_flag.NE.0) GOTO 9000
      nr0b=0; nz0b=0; np0b=0; nfper0=0; nextcur=0
      rminb=0; rmaxb=0; zminb=0; zmaxb=0
      IF (shar_rank.EQ.0) THEN
         CALL cdf_read(ngrid, vn_ir, nr0b, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_ir,
     1      1,REAL((/nr0b/),rprec),ier_flag)
         CALL cdf_read(ngrid, vn_jz, nz0b, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_jz,
     1      1,REAL((/nz0b/),rprec),ier_flag)
         CALL cdf_read(ngrid, vn_kp, np0b, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_kp,
     1      1,REAL((/np0b/),rprec),ier_flag)
         CALL cdf_read(ngrid, vn_nfp, nfper0, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_nfp,
     1      1,REAL((/nfper0/),rprec),ier_flag)
         CALL cdf_read(ngrid, vn_nextcur, nextcur, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_nextcur,
     1      1,REAL((/nextcur/),rprec),ier_flag)
         CALL cdf_read(ngrid, vn_rmin, rminb, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_rmin,
     1      1,(/rminb/),ier_flag)
         CALL cdf_read(ngrid, vn_zmin, zminb, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_zmin,
     1      1,(/zminb/),ier_flag)
         CALL cdf_read(ngrid, vn_rmax, rmaxb, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_rmax,
     1      1,(/rmaxb/),ier_flag)
         CALL cdf_read(ngrid, vn_zmax, zmaxb, ier=istat)
         IF (istat.NE.0) ier_flag=9
         IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_zmax,
     1      1,(/zmaxb/),ier_flag)
         CALL cdf_inquire(ngrid, vn_coilcur, dimlens, ier=istat)
         lraw_local = (istat.EQ.0)
         mgrid_mode = 'N'
         CALL cdf_inquire(ngrid, vn_mgmode, dimlens, ier=istat)
         IF (istat.EQ.0) THEN
            CALL cdf_read(ngrid, vn_mgmode, mgrid_mode, ier=istat)
            IF (istat.NE.0) ier_flag=9
         END IF
      END IF
      CALL mgrid_sync_error(ier_flag, comm)
      IF (ier_flag.NE.0) GOTO 9000
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         CALL MPI_BCAST(nr0b,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(nz0b,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(np0b,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(nfper0,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(nextcur,1,MPI_INTEGER,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(rminb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(zminb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(rmaxb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(zmaxb,1,MPI_REAL8,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(lraw_local,1,MPI_LOGICAL,0,shar_comm,MPI_ERR)
         CALL MPI_BCAST(mgrid_mode,1,MPI_CHARACTER,0,shar_comm,MPI_ERR)
      END IF
#endif
      CALL mgrid_check_header(nv,nfp,SIZE(extcur),ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      delrb=(rmaxb-rminb)/(nr0b-1)
      delzb=(zmaxb-zminb)/(nz0b-1)
      nbvac=nr0b*nz0b*nv

      istat=0
      IF (ASSOCIATED(curlabel)) THEN
         IF (SIZE(curlabel).NE.nextcur) DEALLOCATE(curlabel)
      END IF
      IF (.NOT.ASSOCIATED(curlabel))
     1   ALLOCATE(curlabel(nextcur),STAT=istat)
      IF (istat.NE.0) ier_flag=9
      IF (ASSOCIATED(raw_coil_current)) DEALLOCATE(raw_coil_current)
      IF (lraw_local) THEN
         ALLOCATE(raw_coil_current(nextcur),STAT=istat)
         IF (istat.NE.0) ier_flag=9
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      curlabel=' '
      IF (shar_rank.EQ.0) THEN
         IF (nextcur.EQ.1) THEN
            CALL cdf_read(ngrid,vn_coilgrp,curlabel(1),ier=istat)
         ELSE
            CALL cdf_read(ngrid,vn_coilgrp,curlabel,ier=istat)
         END IF
         IF (istat.NE.0) ier_flag=9
         IF (lraw_local) THEN
            CALL cdf_read(ngrid,vn_coilcur,raw_coil_current,ier=istat)
            IF (istat.NE.0) ier_flag=9
            IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,vn_coilcur,
     1         SIZE(raw_coil_current),raw_coil_current,ier_flag)
         END IF
         CALL cdf_close(ngrid,istat)
         file_open=.FALSE.
         IF (istat.NE.0) ier_flag=9
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         CALL MPI_BCAST(curlabel,SIZE(curlabel)*LEN(curlabel),
     1                  MPI_CHARACTER,0,shar_comm,MPI_ERR)
         IF (lraw_local) CALL MPI_BCAST(raw_coil_current,nextcur,
     1                                MPI_REAL8,0,shar_comm,MPI_ERR)
      END IF
#endif
      CALL mgrid_prepare_bvac(nbvac,shar_rank,shar_comm,
     1                        lmpi_local,ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      IF (shar_rank.EQ.0) bvac=0
#if defined(MPI_OPT)
      IF (lmpi_local) CALL MPI_BARRIER(shar_comm,MPI_ERR)
#endif
!     Preserve the original parallel decomposition over current groups.
      temp_rank=-1
      IF (mpi_rank.LT.nextcur) THEN
         CALL cdf_open(ngrid,filename,'r',istat)
         file_open=(istat.EQ.0)
         IF (istat.NE.0) ier_flag=9
         ALLOCATE(brtemp(nr0b,nz0b,np0b),bptemp(nr0b,nz0b,np0b),
     1            bztemp(nr0b,nz0b,np0b),bttemp(nbvac,3),STAT=istat)
         IF (istat.NE.0) ier_flag=9
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         ig=MPI_UNDEFINED
         IF (mpi_rank.LT.nextcur) ig=0
         CALL MPI_COMM_SPLIT(comm,ig,mpi_rank,temp_comm,istat)
         IF (ig.EQ.0) THEN
            CALL MPI_COMM_RANK(temp_comm,temp_rank,istat)
            CALL MPI_COMM_SIZE(temp_comm,temp_size,istat)
         END IF
      ELSE
#endif
         temp_rank=0
         temp_size=1
#if defined(MPI_OPT)
      END IF
#endif
      IF (temp_rank.GE.0) THEN
         bttemp=0
         DO ig=temp_rank+1,nextcur,temp_size
            WRITE(temp,1000) vn_br0,ig
            CALL cdf_inquire(ngrid,temp,dimlens,ier=istat)
            IF (istat.NE.0) THEN
               ier_flag=9
               EXIT
            END IF
            IF (ANY(dimlens.NE.(/nr0b,nz0b,np0b/))) THEN
               ier_flag=9
               EXIT
            END IF
            CALL cdf_read(ngrid,temp,brtemp,ier=istat)
            IF (istat.NE.0) ier_flag=9
            IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,temp,
     1         SIZE(brtemp),brtemp,ier_flag)
            WRITE(temp,1000) vn_bp0,ig
            CALL cdf_inquire(ngrid,temp,dimlens,ier=istat)
            IF (istat.NE.0) THEN
               ier_flag=9
               EXIT
            END IF
            IF (ANY(dimlens.NE.(/nr0b,nz0b,np0b/))) THEN
               ier_flag=9
               EXIT
            END IF
            CALL cdf_read(ngrid,temp,bptemp,ier=istat)
            IF (istat.NE.0) ier_flag=9
            IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,temp,
     1         SIZE(bptemp),bptemp,ier_flag)
            WRITE(temp,1000) vn_bz0,ig
            CALL cdf_inquire(ngrid,temp,dimlens,ier=istat)
            IF (istat.NE.0) THEN
               ier_flag=9
               EXIT
            END IF
            IF (ANY(dimlens.NE.(/nr0b,nz0b,np0b/))) THEN
               ier_flag=9
               EXIT
            END IF
            CALL cdf_read(ngrid,temp,bztemp,ier=istat)
            IF (istat.NE.0) ier_flag=9
            IF (istat.EQ.0) CALL mgrid_nc_values(ngrid,temp,
     1         SIZE(bztemp),bztemp,ier_flag)
            IF (ier_flag.NE.0) EXIT
            IF (.NOT.ALL(ieee_is_finite(brtemp)) .OR.
     1          .NOT.ALL(ieee_is_finite(bptemp)) .OR.
     2          .NOT.ALL(ieee_is_finite(bztemp))) THEN
               ier_flag=9
               EXIT
            END IF
            CALL sum_bfield(bttemp(:,1),brtemp,extcur(ig),nv)
            CALL sum_bfield(bttemp(:,2),bptemp,extcur(ig),nv)
            CALL sum_bfield(bttemp(:,3),bztemp,extcur(ig),nv)
         END DO
         IF (.NOT.ALL(ieee_is_finite(bttemp))) ier_flag=9
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      IF (ier_flag.NE.0) GOTO 9000
      np0b=nv
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         IF (temp_rank.GE.0) THEN
            CALL MPI_REDUCE(bttemp,bvac,SIZE(bttemp),MPI_REAL8,
     1                      MPI_SUM,0,temp_comm,istat)
            CALL MPI_COMM_FREE(temp_comm,istat)
         END IF
         ig=MPI_UNDEFINED
         IF (shar_rank.EQ.0) ig=0
         CALL MPI_COMM_SPLIT(comm,ig,mpi_rank,temp_comm,istat)
         IF (shar_rank.EQ.0) THEN
            CALL MPI_BCAST(bvac,SIZE(bvac),MPI_REAL8,0,temp_comm,istat)
            CALL MPI_COMM_FREE(temp_comm,istat)
         END IF
      ELSE
#endif
         bvac=bttemp
#if defined(MPI_OPT)
      END IF
#endif
      IF (shar_rank.EQ.0) THEN
         IF (.NOT.ALL(ieee_is_finite(bvac))) ier_flag=9
      END IF
      CALL mgrid_sync_error(ier_flag,comm)
      nobser=0; nobd=0; nbsets=0
!     Every exit closes files and communicators, including partial reads.
 9000 IF (file_open) CALL cdf_close(ngrid,istat)
      IF (ALLOCATED(brtemp)) DEALLOCATE(brtemp)
      IF (ALLOCATED(bptemp)) DEALLOCATE(bptemp)
      IF (ALLOCATED(bztemp)) DEALLOCATE(bztemp)
      IF (ALLOCATED(bttemp)) DEALLOCATE(bttemp)
#if defined(MPI_OPT)
      IF (lmpi_local) THEN
         IF (temp_comm.NE.MPI_COMM_NULL)
     1      CALL MPI_COMM_FREE(temp_comm,istat)
         CALL MPI_COMM_FREE(shar_comm,istat)
      END IF
#endif
      IF (ier_flag.EQ.0) CALL assign_bptrs(bvac)
 1000 FORMAT(a,'_',i3.3)
      END SUBROUTINE read_mgrid_nc

!     NetCDF missing-data markers are often finite. Never interpret
!     unwritten field samples or declared missing values as coil data.
!     Check exact markers, without imposing an artificial magnitude cut.
      SUBROUTINE mgrid_nc_values(ncid,name,nvalues,values,ier)
      USE netcdf, ONLY: nf90_inq_varid, nf90_inquire_variable,
     1   nf90_inquire_attribute, nf90_get_att, nf90_noerr, nf90_enotatt,
     2   nf90_double, nf90_int, nf90_fill_double, nf90_fill_int
      INTEGER, INTENT(in) :: ncid,nvalues
      CHARACTER(LEN=*), INTENT(in) :: name
      REAL(rprec), INTENT(in) :: values(nvalues)
      INTEGER, INTENT(inout) :: ier
      INTEGER :: varid,xtype,atttype,nmarkers,status,i,j
      REAL(rprec) :: default_fill
      REAL(rprec), ALLOCATABLE :: markers(:)
      CHARACTER(LEN=13), PARAMETER :: attributes(2) =
     1   [CHARACTER(LEN=13) :: '_FillValue','missing_value']

      IF (ier.NE.0) RETURN
      IF (.NOT.ALL(ieee_is_finite(values))) THEN
         ier=9
         RETURN
      END IF
      status=nf90_inq_varid(ncid,TRIM(name),varid)
      IF (status.NE.nf90_noerr) GOTO 9000
      status=nf90_inquire_variable(ncid,varid,xtype=xtype)
      IF (status.NE.nf90_noerr) GOTO 9000
!     These are the existing EZCDF types used by MGRID numeric reads.
      SELECT CASE(xtype)
      CASE(nf90_double)
         default_fill=nf90_fill_double
      CASE(nf90_int)
         default_fill=REAL(nf90_fill_int,rprec)
      CASE DEFAULT
         GOTO 9000
      END SELECT
      DO i=1,SIZE(attributes)
         status=nf90_inquire_attribute(ncid,varid,
     1      TRIM(attributes(i)),xtype=atttype,len=nmarkers)
         IF (status.EQ.nf90_enotatt) THEN
            IF (i.EQ.1) THEN
               IF (ANY(values.EQ.default_fill)) GOTO 9000
            END IF
            CYCLE
         END IF
         IF (status.NE.nf90_noerr) GOTO 9000
         IF (nmarkers.LT.1) GOTO 9000
!        _FillValue is a scalar with the variable's storage type.
         IF (i.EQ.1) THEN
            IF (nmarkers.NE.1 .OR. atttype.NE.xtype) GOTO 9000
         END IF
         ALLOCATE(markers(nmarkers),STAT=status)
         IF (status.NE.0) GOTO 9000
         status=nf90_get_att(ncid,varid,TRIM(attributes(i)),markers)
         IF (status.NE.nf90_noerr) GOTO 9000
         DO j=1,nmarkers
!           Nonfinite samples were rejected above, including NaN fills.
            IF (.NOT.ieee_is_finite(markers(j))) CYCLE
            IF (ANY(values.EQ.markers(j))) GOTO 9000
         END DO
         DEALLOCATE(markers)
      END DO
      RETURN
 9000 ier=9
      END SUBROUTINE mgrid_nc_values
#endif

      SUBROUTINE sum_bfield(bfield, bf_add, cur, nv)
      INTEGER, INTENT(IN)        :: nv
      REAL(rprec), INTENT(INOUT) :: bfield(nr0b*nz0b,nv)
      REAL(rprec), INTENT(IN)    :: bf_add(nr0b*nz0b,np0b)
      INTEGER     :: nskip
      REAL(rprec) :: cur
      nskip = np0b/nv
      bfield = bfield + cur*bf_add(:,1:np0b:nskip)

      END SUBROUTINE sum_bfield

      SUBROUTINE assign_bptrs(bptr)
      IMPLICIT NONE
      REAL(rprec), TARGET, INTENT(in) :: bptr(nr0b,nz0b,np0b,3)

      brvac => bptr(:,:,:,1)
      bpvac => bptr(:,:,:,2)
      bzvac => bptr(:,:,:,3)

      END SUBROUTINE assign_bptrs

      SUBROUTINE free_mgrid_diagnostics
!     Optional sections can be absent or only partly allocated.
      IF (ASSOCIATED(needflx)) DEALLOCATE(needflx)
      IF (ASSOCIATED(nbcoils)) DEALLOCATE(nbcoils)
      IF (ASSOCIATED(limitr)) DEALLOCATE(limitr)
      IF (ASSOCIATED(nsetsn)) DEALLOCATE(nsetsn)
      IF (ASSOCIATED(iconnect)) DEALLOCATE(iconnect)
      IF (ASSOCIATED(needbfld)) DEALLOCATE(needbfld)
      IF (ASSOCIATED(unpsiext)) DEALLOCATE(unpsiext)
      IF (ASSOCIATED(plbfld)) DEALLOCATE(plbfld)
      IF (ASSOCIATED(rbcoil)) DEALLOCATE(rbcoil)
      IF (ASSOCIATED(zbcoil)) DEALLOCATE(zbcoil)
      IF (ASSOCIATED(abcoil)) DEALLOCATE(abcoil)
      IF (ASSOCIATED(bcoil)) DEALLOCATE(bcoil)
      IF (ASSOCIATED(rbcoilsqr)) DEALLOCATE(rbcoilsqr)
      IF (ASSOCIATED(xobser)) DEALLOCATE(xobser)
      IF (ASSOCIATED(zobser)) DEALLOCATE(zobser)
      IF (ASSOCIATED(xobsqr)) DEALLOCATE(xobsqr)
      IF (ASSOCIATED(dsiext)) DEALLOCATE(dsiext)
      IF (ASSOCIATED(psiext)) DEALLOCATE(psiext)
      IF (ASSOCIATED(plflux)) DEALLOCATE(plflux)
      IF (ASSOCIATED(b_chi)) DEALLOCATE(b_chi)
      IF (ASSOCIATED(dsilabel)) DEALLOCATE(dsilabel)
      IF (ASSOCIATED(bloopnames)) DEALLOCATE(bloopnames)
      IF (ASSOCIATED(dbcoil)) DEALLOCATE(dbcoil)
      IF (ASSOCIATED(pfcspec)) DEALLOCATE(pfcspec)
      IF (ASSOCIATED(rlim)) DEALLOCATE(rlim)
      IF (ASSOCIATED(zlim)) DEALLOCATE(zlim)
      IF (ASSOCIATED(reslim)) DEALLOCATE(reslim)
      IF (ASSOCIATED(seplim)) DEALLOCATE(seplim)
      nobser=0; nobd=0; nbsets=0; nbcoil_max=0; nbfldn=0
      nbcoilsn=0; nlim=0; nlim_max=0; nsets=0
      nrgrid=0; nzgrid=0; tokid=' '
      rx1=0; rx2=0; zy1=0; zy2=0; condif=0
      END SUBROUTINE free_mgrid_diagnostics

      SUBROUTINE free_mgrid (istat,comm)
      USE mpi_inc
      INTEGER, INTENT(out) :: istat
      INTEGER, INTENT(in), OPTIONAL :: comm
#if defined(MPI_OPT)
      LOGICAL :: initialized, collective, same_group, complete
      INTEGER :: ierr
#endif

      istat = 0
#if defined(MPI_OPT)
      collective = .FALSE.
      IF (PRESENT(comm)) THEN
         CALL MPI_INITIALIZED(initialized,ierr)
         collective = initialized .AND. ierr.EQ.MPI_SUCCESS
      END IF
      IF (collective) THEN
         CALL mgrid_window_group(comm,same_group,complete,istat)
         IF (.NOT.complete) istat = 9
         CALL mgrid_sync_error(istat,comm)
         IF (istat.NE.0) RETURN
         CALL mgrid_release_shared_bvac(istat)
         CALL mgrid_sync_error(istat,comm)
         IF (istat.NE.0) RETURN
      ELSE IF (bvac_is_shared) THEN
!        Without a communicator the caller cannot free a live window.
         istat = 9
         RETURN
      END IF
#endif
      IF (ASSOCIATED(bvac)) DEALLOCATE (bvac,stat=istat)
      NULLIFY(brvac,bpvac,bzvac)
      CALL free_mgrid_diagnostics
      IF (ASSOCIATED(curlabel)) DEALLOCATE(curlabel)
      IF (ASSOCIATED(raw_coil_current)) DEALLOCATE(raw_coil_current)

      CALL mgrid_invalidate_cache
      mgrid_nv_old = -1
      mgrid_nfp_old = -1
      mgrid_nextcur_old = 0
      mgrid_mode_old = ' '
      mgrid_identity_old = ' '
      mgrid_stamp_old = 0
      IF (ALLOCATED(mgrid_extcur_old)) DEALLOCATE(mgrid_extcur_old)

      END SUBROUTINE free_mgrid

      END MODULE mgrid_mod
