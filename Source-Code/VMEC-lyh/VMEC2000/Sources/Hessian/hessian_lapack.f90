module hessian_lapack
   use stel_kinds, only: dp
#if defined(VMEC_OPENBLAS)
   use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_size_t, &
                                       c_int, c_null_char, c_associated
#endif
   implicit none
   private
   public :: hessian_getrf, hessian_getrs, use_unblocked_hessian_lu

#if defined(VMEC_OPENBLAS)
   interface
      function openblas_get_config() bind(C) result(value)
         import :: c_ptr
         type(c_ptr) :: value
      end function

      function openblas_get_corename() bind(C) result(value)
         import :: c_ptr
         type(c_ptr) :: value
      end function

      function strncmp(lhs, rhs, n) bind(C) result(order)
         import :: c_ptr, c_char, c_size_t, c_int
         type(c_ptr), value :: lhs
         character(kind=c_char), intent(in) :: rhs(*)
         integer(c_size_t), value :: n
         integer(c_int) :: order
      end function
   end interface
#endif

contains

   logical function use_unblocked_hessian_lu()
      ! The affected OpenBLAS dispatch can corrupt blocked LU and
      ! multi-RHS solves while returning INFO=0. Keep other backends.
#if defined(VMEC_OPENBLAS)
      character(kind=c_char,len=*), parameter :: version = 'OpenBLAS 0.3.20 '
      character(kind=c_char,len=*), parameter :: core = 'Cooperlake'//c_null_char
      logical, save :: checked = .false., affected = .false.
      type(c_ptr) :: config, corename

      if (.not.checked) then
         config = openblas_get_config()
         corename = openblas_get_corename()
         if (c_associated(config) .and. c_associated(corename)) then
            affected = strncmp(config, version, len(version,kind=c_size_t)) == 0 .and. &
                       strncmp(corename, core, len(core,kind=c_size_t)) == 0
         endif
         checked = .true.
      endif
      use_unblocked_hessian_lu = affected
#else
      use_unblocked_hessian_lu = .false.
#endif
   end function

   subroutine hessian_getrf(n, a, piv, info)
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n,n)
      integer, intent(out) :: piv(n), info

      if (use_unblocked_hessian_lu()) then
         call dgetf2(n, n, a, n, piv, info)
      else
         call dgetrf(n, n, a, n, piv, info)
      endif
   end subroutine

   subroutine hessian_getrs(n, nrhs, a, piv, b, info)
      integer, intent(in) :: n, nrhs, piv(n)
      real(dp), intent(in) :: a(n,n)
      real(dp), intent(inout) :: b(n,nrhs)
      integer, intent(out) :: info
      integer :: column

      if (nrhs > 1 .and. use_unblocked_hessian_lu()) then
         do column = 1, nrhs
            call dgetrs('N', n, 1, a, n, piv, b(:,column), n, info)
            if (info /= 0) return
         enddo
      else
         call dgetrs('N', n, nrhs, a, n, piv, b, n, info)
      endif
   end subroutine

end module hessian_lapack
