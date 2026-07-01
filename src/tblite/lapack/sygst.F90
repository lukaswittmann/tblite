! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/lapack/sygst.F90
!> Reduces a real symmetric-definite generalized eigenproblem to standard form.

! Integer kind used to interface with the external BLAS/LAPACK backend.
! Defaults to 32-bit indices (LP64); define IK=i8 to build against an
! ILP64 (64-bit integer) linear algebra backend.
#ifndef IK
#define IK i4
#endif

!> Reduces a real symmetric-definite generalized eigenproblem to standard form.
module tblite_lapack_sygst
   use mctc_env, only : sp, dp, ik => IK
   implicit none
   private

   public :: wrap_sygst


   !> Reduces a real symmetric-definite generalized eigenproblem to standard form.
   !>
   !> If ITYPE = 1, the problem is A*x = lambda*B*x,
   !> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !>
   !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !>
   !> B must have been previously factorized as U**T*U or L*L**T by POTRF.
   interface wrap_sygst
      module procedure :: wrap_ssygst
      module procedure :: wrap_dsygst
   end interface wrap_sygst


   !> Reduces a real symmetric-definite generalized eigenproblem to standard form.
   !>
   !> If ITYPE = 1, the problem is A*x = lambda*B*x,
   !> and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
   !>
   !> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
   !> B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
   !>
   !> B must have been previously factorized as U**T*U or L*L**T by POTRF.
   interface lapack_sygst
      pure subroutine ssygst(itype, uplo, n, a, lda, b, ldb, info)
         import :: sp, ik
         integer(ik), intent(in) :: itype
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
      end subroutine ssygst
      pure subroutine dsygst(itype, uplo, n, a, lda, b, ldb, info)
         import :: dp, ik
         integer(ik), intent(in) :: itype
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         character(len=1), intent(in) :: uplo
      end subroutine dsygst
   end interface lapack_sygst

contains

pure subroutine wrap_ssygst(amat, bmat, info, itype, uplo)
   real(sp), intent(inout) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   integer, intent(out) :: info
   character(len=1) :: ula
   integer(ik) :: ita, n, lda, ldb, stat

   ita = 1
   if(present(itype)) ita = itype
   ula = 'u'
   if(present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   call lapack_sygst(ita, ula, n, amat, lda, bmat, ldb, stat)
   info = stat
end subroutine wrap_ssygst

pure subroutine wrap_dsygst(amat, bmat, info, itype, uplo)
   real(dp), intent(inout) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   integer, intent(in), optional :: itype
   character(len=1), intent(in), optional :: uplo
   integer, intent(out) :: info
   character(len=1) :: ula
   integer(ik) :: ita, n, lda, ldb, stat

   ita = 1
   if(present(itype)) ita = itype
   ula = 'u'
   if(present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   call lapack_sygst(ita, ula, n, amat, lda, bmat, ldb, stat)
   info = stat
end subroutine wrap_dsygst

end module tblite_lapack_sygst
