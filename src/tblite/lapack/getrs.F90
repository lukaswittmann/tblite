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

!> @file tblite/lapack/getrs.F90
!> Provides wrappers to solve a linear equation system

! Integer kind used to interface with the external BLAS/LAPACK backend.
! Defaults to 32-bit indices (LP64); define IK=i8 to build against an
! ILP64 (64-bit integer) linear algebra backend.
#ifndef IK
#define IK i4
#endif

!> Wrappers to solve a system of linear equations
module tblite_lapack_getrs
   use mctc_env, only : sp, dp, ik => IK
   implicit none
   private

   public :: wrap_getrs


   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface wrap_getrs
      module procedure :: wrap_sgetrs
      module procedure :: wrap_dgetrs
   end interface wrap_getrs


   !> Solves a system of linear equations
   !>    A * X = B  or  A**T * X = B
   !> with a general N-by-N matrix A using the LU factorization computed
   !> by ?GETRF.
   interface lapack_getrs
      pure subroutine sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp, ik
         integer(ik), intent(in) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
      end subroutine sgetrs
      pure subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp, ik
         integer(ik), intent(in) :: ipiv(*)
         integer(ik), intent(out) :: info
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: nrhs
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         character(len=1), intent(in) :: trans
      end subroutine dgetrs
   end interface lapack_getrs

contains

subroutine wrap_sgetrs(amat, bmat, ipiv, info, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer(ik) :: n, nrhs, lda, ldb, stat
   integer(ik) :: jpiv(size(ipiv))
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   jpiv(:) = ipiv
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, jpiv, bmat, ldb, stat)
   info = stat
end subroutine wrap_sgetrs


subroutine wrap_dgetrs(amat, bmat, ipiv, info, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   integer, intent(out) :: info
   character(len=1), intent(in), optional :: trans
   character(len=1) :: tra
   integer(ik) :: n, nrhs, lda, ldb, stat
   integer(ik) :: jpiv(size(ipiv))
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   jpiv(:) = ipiv
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_getrs(tra, n, nrhs, amat, lda, jpiv, bmat, ldb, stat)
   info = stat
end subroutine wrap_dgetrs

end module tblite_lapack_getrs
