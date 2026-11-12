      MODULE vmec_io
      USE stel_kinds
      USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_char,
     &                                     c_null_char
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: VolAvgB, IonLarmor, Aminor_p, Rmajor_p,
     1  betatot, betapol, betator, betaxis, b0, volume_p,
     2  cross_area_p, surf_area_p, circum_p, kappa_p
      REAL(rprec) :: rmax_surf, rmin_surf, zmax_surf
      PRIVATE :: c_int, c_char, c_null_char, c_mkstemps,
     &           c_close, c_rename, c_unlink

      INTERFACE
         INTEGER(c_int) FUNCTION c_mkstemps(path, suffix_len)
     &      BIND(C,NAME='mkstemps')
         IMPORT c_int, c_char
         CHARACTER(KIND=c_char), INTENT(INOUT) :: path(*)
         INTEGER(c_int), VALUE :: suffix_len
         END FUNCTION

         INTEGER(c_int) FUNCTION c_close(fd) BIND(C,NAME='close')
         IMPORT c_int
         INTEGER(c_int), VALUE :: fd
         END FUNCTION

         INTEGER(c_int) FUNCTION c_rename(oldpath, newpath)
     &      BIND(C,NAME='rename')
         IMPORT c_int, c_char
         CHARACTER(KIND=c_char), INTENT(IN) :: oldpath(*), newpath(*)
         END FUNCTION

         INTEGER(c_int) FUNCTION c_unlink(path) BIND(C,NAME='unlink')
         IMPORT c_int, c_char
         CHARACTER(KIND=c_char), INTENT(IN) :: path(*)
         END FUNCTION
      END INTERFACE

      CONTAINS

      SUBROUTINE reserve_output_path(filename, tag, path, status)
      CHARACTER(LEN=*), INTENT(IN) :: filename, tag
      CHARACTER(LEN=*), INTENT(OUT) :: path
      INTEGER, INTENT(OUT) :: status
      CHARACTER(KIND=c_char,LEN=:), ALLOCATABLE :: name
      INTEGER :: dot
      INTEGER(c_int) :: fd

      path = ''
      status = 1
      dot = INDEX(filename, '.', BACK=.TRUE.)
      IF (dot.EQ.0) dot = LEN_TRIM(filename)+1
      name = filename(:dot-1)//tag//'XXXXXX'//
     &       TRIM(filename(dot:))//c_null_char
      IF (LEN(name)-1.GT.LEN(path)) RETURN

      fd = c_mkstemps(name, INT(LEN_TRIM(filename)-dot+1,c_int))
      IF (fd.LT.0) RETURN
      path = name(:LEN(name)-1)
      status = c_close(fd)
      IF (status.NE.0) CALL discard_output_path(path)
      END SUBROUTINE reserve_output_path

      SUBROUTINE publish_output_path(path, filename, status)
      CHARACTER(LEN=*), INTENT(IN) :: path, filename
      INTEGER, INTENT(OUT) :: status

      status = c_rename(TRIM(path)//c_null_char,
     &                  TRIM(filename)//c_null_char)
      IF (status.NE.0) PRINT *, 'Cannot publish ', TRIM(filename)
      END SUBROUTINE publish_output_path

      SUBROUTINE discard_output_path(path)
      CHARACTER(LEN=*), INTENT(INOUT) :: path
      INTEGER(c_int) :: status

      IF (LEN_TRIM(path).EQ.0) RETURN
      status = c_unlink(TRIM(path)//c_null_char)
      path = ''
      END SUBROUTINE discard_output_path

      SUBROUTINE equilibrium_output_names(extension, names)
      CHARACTER(LEN=*), INTENT(IN) :: extension
      CHARACTER(LEN=*), INTENT(OUT) :: names(3)

#ifdef NETCDF
      names(1) = 'jxbout_'//TRIM(extension)//'.nc'
#else
      names(1) = 'jxbout.'//TRIM(extension)//'.txt'
#endif
      names(2) = 'wout_'//TRIM(extension)//'.txt'
      names(3) = 'wout_'//TRIM(extension)//'.nc'
      END SUBROUTINE equilibrium_output_names

      SUBROUTINE commit_equilibrium_outputs(extension, paths, status)
      CHARACTER(LEN=*), INTENT(IN) :: extension
      CHARACTER(LEN=*), INTENT(INOUT) :: paths(3)
      INTEGER, INTENT(OUT) :: status
      CHARACTER(LEN=128) :: names(3)
      LOGICAL :: published(3)
      INTEGER :: i

      CALL equilibrium_output_names(extension,names)
      published = .FALSE.
      status = 0
!     Publish the canonical NetCDF WOUT only after its diagnostics.
      DO i = 1, SIZE(paths)
         IF (LEN_TRIM(paths(i)).EQ.0) CYCLE
         CALL publish_output_path(paths(i), names(i), status)
         IF (status.NE.0) EXIT
         published(i) = .TRUE.
         paths(i) = ''
      END DO
      IF (status.NE.0) THEN
         DO i = 1, SIZE(paths)
            IF (published(i)) CALL discard_output_path(names(i))
         END DO
      END IF
      END SUBROUTINE commit_equilibrium_outputs

      SUBROUTINE archive_equilibrium_outputs(extension, reset, status)
      CHARACTER(LEN=*), INTENT(IN) :: extension
      CHARACTER(LEN=*), INTENT(INOUT) :: reset
      INTEGER, INTENT(OUT) :: status
      CHARACTER(LEN=128) :: names(3)
      CHARACTER(LEN=160) :: previous
      CHARACTER(LEN=LEN(reset)+10) :: seed
      LOGICAL :: exists, isdir, seed_existed, seed_exists, seed_nc
      INTEGER :: i

      CALL equilibrium_output_names(extension,names)
      status = 0
      seed = 'wout'
      seed_existed = .FALSE.
      IF (LEN_TRIM(reset).GT.0) THEN
         CALL parse_extension(seed, reset, seed_nc)
         INQUIRE(FILE=TRIM(seed), EXIST=seed_existed)
      END IF
      DO i = 1, SIZE(names)
         INQUIRE(FILE=TRIM(names(i)), EXIST=exists, IOSTAT=status)
         IF (status.NE.0) RETURN
         IF (.NOT.exists) CYCLE
         INQUIRE(FILE=TRIM(names(i))//'/.', EXIST=isdir)
         IF (isdir) CYCLE
         CALL reserve_output_path(names(i),'.previous.',previous,
     &                            status)
         IF (status.NE.0) RETURN
         CALL publish_output_path(names(i), previous, status)
         IF (status.NE.0) THEN
            CALL discard_output_path(previous)
            RETURN
         END IF
         PRINT *, 'Previous equilibrium output saved as ',
     &            TRIM(previous)
         IF (i.GE.2 .AND. seed_existed) THEN
            INQUIRE(FILE=TRIM(seed), EXIST=seed_exists)
            IF (.NOT.seed_exists) THEN
               reset = previous
               seed = previous
            END IF
         END IF
      END DO
      END SUBROUTINE archive_equilibrium_outputs
      END MODULE vmec_io
