      SUBROUTINE becoil (rad, zee, br, bp, bz, brvac, bpvac, bzvac,            &
     &                   ier_flag)
      USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
      USE vmec_params, ONLY: norm_term_flag, vacuum_error_flag
      USE vacmod
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nuv3), INTENT(in) :: rad, zee
      REAL(dp), DIMENSION(nuv3), INTENT(out) :: br, bp, bz
      REAL(dp), DIMENSION(nr0b,nz0b,np0b), INTENT(in) ::
     1   brvac, bpvac, bzvac
      INTEGER, INTENT(out) :: ier_flag
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, kv, ir, jz, ir1, jz1
      REAL(dp) :: rad0, zee0, ri, zj,
     1   pr, qz, w22, w21, w12, w11, tbecon, tbecoff
C-----------------------------------------------
!
!     DETERMINE THE CYLINDRICAL COMPONENTS OF THE EXTERNAL
!     MAGNETIC FIELD (BR, BP, BZ) AT A FIXED PHI PLANE BY 
!     USING 2-D INTERPOLATION BASED ON THE FOUR POINT FORMULA
!     IN ABRAMOWITZ AND STEGUN, EQ. 25.2.66 
!
!     BRVAC, BPVAC, BZVAC:       CYLINDRICAL COMPONENTS OF VACUUM B-FIELD 
!                                STORED ON R, Z, PHI GRID
!
!     RAD, ZEE:                  ARRAY OF R, Z VALUES (ON FLUX SURFACE)
!                                AT WHICH B-FIELD IS TO BE DETERMINED
!
C-----------------------------------------------
      CALL second0(tbecon)

!     SURFACE reconstructs all NESTOR sample coordinates on each rank.
!     Reject an uncovered boundary instead of substituting grid-edge fields.
      ier_flag = vacuum_error_flag
      IF (.NOT.ALL(ieee_is_finite(rad))) RETURN
      IF (.NOT.ALL(ieee_is_finite(zee))) RETURN
      IF (MINVAL(rad).LT.rminb .OR. MAXVAL(rad).GT.rmaxb) RETURN
      IF (MINVAL(zee).LT.zminb .OR. MAXVAL(zee).GT.zmaxb) RETURN
      ier_flag = norm_term_flag

      DO i = nuv3min, nuv3max
         rad0 = rad(i)
         zee0 = zee(i)
!
!       DETERMINE PHI-PLANE, KV (MUST LIE IN FIRST FIELD PERIOD)
!
         kv = 1 + MOD(i - 1,nv)
         kv = MIN (kv, np0b)                       !!Axi-symmetric special case
!
!
!       DETERMINE INTEGER INDICES (IR,JZ) FOR LOWER LEFT R, Z CORNER GRID POINT
!
         ir = MIN(nr0b, INT((rad0 - rminb)/delrb) + 1)
         jz = MIN(nz0b, INT((zee0 - zminb)/delzb) + 1)
         ir1 = MIN(nr0b,ir + 1)
         jz1 = MIN(nz0b,jz + 1)
!
!       COMPUTE RI, ZJ AND PR , QZ AT GRID POINT (IR , JZ)
!       ALSO, COMPUTE WEIGHTS WIJ FOR 4 CORNER GRID POINTS
!
         ri = rminb + (ir - 1)*delrb
         zj = zminb + (jz - 1)*delzb
         pr = (rad0 - ri)/delrb
         qz = (zee0 - zj)/delzb
         w22 = pr*qz                   !p*q               
         w21 = pr - w22                !p*(1-q)
         w12 = qz - w22                !q*(1-p)
         w11 = 1 + w22 - (pr + qz)     !(1-p)*(1-q)
!
!       COMPUTE B FIELD AT R, PHI, Z BY INTERPOLATION
!
         br(i) = w11*brvac(ir,jz,kv) + w22*brvac(ir1,jz1,kv) +
     1      w21*brvac(ir1,jz,kv) + w12*brvac(ir,jz1,kv)
         bz(i) = w11*bzvac(ir,jz,kv) + w22*bzvac(ir1,jz1,kv) +
     1      w21*bzvac(ir1,jz,kv) + w12*bzvac(ir,jz1,kv)
         bp(i) = w11*bpvac(ir,jz,kv) + w22*bpvac(ir1,jz1,kv) +
     1      w21*bpvac(ir1,jz,kv) + w12*bpvac(ir,jz1,kv)

      END DO

      CALL second0(tbecoff)
      becoil_time = becoil_time + (tbecoff - tbecon)

      END SUBROUTINE becoil
