//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef KOKKOSBLAS3_TRSM_HPP_
#define KOKKOSBLAS3_TRSM_HPP_

/// \file KokkosBlas3_trsm.hpp

#include "KokkosKernels_Macros.hpp"
#include "KokkosBlas3_trsm_spec.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Error.hpp"
#include <sstream>
#include <type_traits>

namespace KokkosBlas {

/// \brief Solve triangular linear system with multiple RHSs:
///        op(A)*X = alpha*B if side == "L" or "l"
///        X*op(A) = alpha*B if side == "R" or "r"
/// This function is currently blocking when running the native implementation
/// which only has a serial implementation.
///
/// \tparam execution_space a Kokkos execution space to run the kernels on.
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
/// \tparam BViewType Input(RHS)/Output(solution) M-by-N matrix, as a 2-D
/// Kokkos::View
///
/// \param space [in] an execution space instance that may contain a stream
/// or a queue to execute the kernel on, this only works with TPLs at the
/// moment. \param side  [in] "L" or "l" indicates matrix A is on the left of X
///                   "R" or "r" indicates matrix A is on the right of X
/// \param uplo  [in] "U" or "u" indicates matrix A upper part is stored, the
/// other part is not referenced
///                   "L" or "l" indicates matrix A lower part is stored, the
///                   other part is not referenced
/// \param trans [in] "N" or "n" for non-transpose, "T" or "t" for transpose,
/// "C" or "c" for conjugate transpose.
/// \param diag  [in] "U" or "u" indicates the diagonal of A is assumed to be
/// unit
//                    "N" or "n" indicated the diagonal of A is assumed to be
//                    non-unit
/// \param alpha [in] Input coefficient used for multiplication with B
/// \param A [in]     Input matrix, as a 2-D Kokkos::View
///                   If side == "L" or "l", matrix A is a M-by-M triangular
///                   matrix; otherwise, matrix A is a N-by-N triangular matrix
/// \param B [in,out] Input/Output matrix, as a 2-D Kokkos::View
///                   On entry, M-by-N matrix of multile RHS
///                   On exit, overwritten with the solution X
template <class execution_space, class AViewType, class BViewType>
void trsm(const execution_space& space, const char side[], const char uplo[], const char trans[], const char diag[],
          typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B) {
  static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<BViewType>::value, "BViewType must be a Kokkos::View.");
  static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");
  static_assert(static_cast<int>(BViewType::rank) == 2, "BViewType must have rank 2.");

  // Check validity of indicator argument
  bool valid_side  = (side[0] == 'L') || (side[0] == 'l') || (side[0] == 'R') || (side[0] == 'r');
  bool valid_uplo  = (uplo[0] == 'U') || (uplo[0] == 'u') || (uplo[0] == 'L') || (uplo[0] == 'l');
  bool valid_trans = (trans[0] == 'N') || (trans[0] == 'n') || (trans[0] == 'T') || (trans[0] == 't') ||
                     (trans[0] == 'C') || (trans[0] == 'c');
  bool valid_diag = (diag[0] == 'U') || (diag[0] == 'u') || (diag[0] == 'N') || (diag[0] == 'n');
  if (!valid_side) {
    std::ostringstream os;
    os << "KokkosBlas::trsm: side = '" << side[0] << "'. "
       << "Valid values include 'L' or 'l' (A is on the left of X), "
          "'R' or 'r' (A is on the right of X).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (!valid_uplo) {
    std::ostringstream os;
    os << "KokkosBlas::trsm: uplo = '" << uplo[0] << "'. "
       << "Valid values include 'U' or 'u' (A is upper triangular), "
          "'L' or 'l' (A is lower triangular).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (!valid_trans) {
    std::ostringstream os;
    os << "KokkosBlas::trsm: trans = '" << trans[0] << "'. "
       << "Valid values include 'N' or 'n' (No transpose), 'T' or 't' "
          "(Transpose), "
          "and 'C' or 'c' (Conjugate transpose).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (!valid_diag) {
    std::ostringstream os;
    os << "KokkosBlas::trsm: diag = '" << diag[0] << "'. "
       << "Valid values include 'U' or 'u' (the diagonal of A is assumed to be "
          "unit), "
          "'N' or 'n' (the diagonal of A is assumed to be non-unit).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Check compatibility of dimensions at run time.
  bool A_s = (side[0] == 'L' || side[0] == 'l');

  int64_t A0 = A.extent(0);
  int64_t A1 = A.extent(1);
  int64_t B0 = B.extent(0);
  int64_t B1 = B.extent(1);

  if ((A0 != A1) || ((A_s ? B0 : B1) != A1)) {
    std::ostringstream os;
    os << "KokkosBlas::trsm: Dimensions of A and B do not match: "
       << "side: " << side[0] << " A: " << A.extent(0) << " x " << A.extent(1) << " B: " << B.extent(0) << " x "
       << B.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Return if degenerated matrices are provided
  if ((A.extent(0) == 0) || (A.extent(1) == 0) || (B.extent(0) == 0) || (B.extent(1) == 0)) return;

  // Minimize the number of Impl::TRSM instantiations, by
  // standardizing on particular View specializations for its template
  // parameters.
  using AVT = Kokkos::View<typename AViewType::const_value_type**, typename AViewType::array_layout,
                           typename AViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using BVT = Kokkos::View<typename BViewType::non_const_value_type**, typename BViewType::array_layout,
                           typename BViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  KokkosBlas::Impl::TRSM<execution_space, AVT, BVT>::trsm(space, side, uplo, trans, diag, alpha, A, B);
}

/// \brief Solve triangular linear system with multiple RHSs:
///        op(A)*X = alpha*B if side == "L" or "l"
///        X*op(A) = alpha*B if side == "R" or "r"
///
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
/// \tparam BViewType Input(RHS)/Output(solution) M-by-N matrix, as a 2-D
/// Kokkos::View
///
/// \param side  [in] "L" or "l" indicates matrix A is on the left of X
///                   "R" or "r" indicates matrix A is on the right of X
/// \param uplo  [in] "U" or "u" indicates matrix A upper part is stored, the
/// other part is not referenced
///                   "L" or "l" indicates matrix A lower part is stored, the
///                   other part is not referenced
/// \param trans [in] "N" or "n" for non-transpose, "T" or "t" for transpose,
/// "C" or "c" for conjugate transpose.
/// \param diag  [in] "U" or "u" indicates the diagonal of A is assumed to be
/// unit
//                    "N" or "n" indicated the diagonal of A is assumed to be
//                    non-unit
/// \param alpha [in] Input coefficient used for multiplication with B
/// \param A [in]     Input matrix, as a 2-D Kokkos::View
///                   If side == "L" or "l", matrix A is a M-by-M triangular
///                   matrix; otherwise, matrix A is a N-by-N triangular matrix
/// \param B [in,out] Input/Output matrix, as a 2-D Kokkos::View
///                   On entry, M-by-N matrix of multile RHS
///                   On exit, overwritten with the solution X
template <class AViewType, class BViewType>
void trsm(const char side[], const char uplo[], const char trans[], const char diag[],
          typename BViewType::const_value_type& alpha, const AViewType& A, const BViewType& B) {
  trsm(typename AViewType::execution_space{}, side, uplo, trans, diag, alpha, A, B);
}
}  // namespace KokkosBlas

#endif  // KOKKOS_BLAS3_TRSM_HPP_
