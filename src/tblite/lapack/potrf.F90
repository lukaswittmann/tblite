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

!> @file tblite/lapack/potrf.F90
!> Computes the Cholesky factorization of a real symmetric positive definite matrix A.

! Integer kind used to interface with the external BLAS/LAPACK backend.
! Defaults to 32-bit indices (LP64); define IK=i8 to build against an
! ILP64 (64-bit integer) linear algebra backend.
#ifndef IK
#define IK i4
#endif

!> Computes the Cholesky factorization of a real symmetric positive definite matrix A.
module tblite_lapack_potrf
   use mctc_env, only : sp, dp, ik => IK
   implicit none
   private

   public :: wrap_potrf


   !> Computes the Cholesky factorization of a real symmetric
   !> positive definite matrix A.
   !>
   !> The factorization has the form
   !>    A = U**T * U,  if UPLO = 'U', or
   !>    A = L  * L**T,  if UPLO = 'L',
   !> where U is an upper triangular matrix and L is lower triangular.
   !>
   !> This is the block version of the algorithm, calling Level 3 BLAS.
   interface wrap_potrf
      module procedure :: wrap_spotrf
      module procedure :: wrap_dpotrf
   end interface wrap_potrf

   !> Computes the Cholesky factorization of a real symmetric
   !> positive definite matrix A.
   !>
   !> The factorization has the form
   !>    A = U**T * U,  if UPLO = 'U', or
   !>    A = L  * L**T,  if UPLO = 'L',
   !> where U is an upper triangular matrix and L is lower triangular.
   !>
   !> This is the block version of the algorithm, calling Level 3 BLAS.
   interface lapack_potrf
      pure subroutine spotrf(uplo, n, a, lda, info)
         import :: sp, ik
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
      end subroutine spotrf
      pure subroutine dpotrf(uplo, n, a, lda, info)
         import :: dp, ik
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: lda
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
      end subroutine dpotrf
   end interface lapack_potrf

contains

subroutine wrap_spotrf(amat, info, uplo)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: uplo
   integer(ik) :: n, lda, stat
   character(len=1) :: ula

   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, stat)
   info = stat
end subroutine wrap_spotrf

subroutine wrap_dpotrf(amat, info, uplo)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: uplo
   integer(ik) :: n, lda, stat
   character(len=1) :: ula

   ula = 'u'
   if (present(uplo)) ula = uplo
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call lapack_potrf(ula, n, amat, lda, stat)
   info = stat
end subroutine wrap_dpotrf

end module tblite_lapack_potrf
