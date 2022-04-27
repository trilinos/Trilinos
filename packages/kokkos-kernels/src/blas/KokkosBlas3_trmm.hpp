/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSBLAS3_TRMM_HPP_
#define KOKKOSBLAS3_TRMM_HPP_

/// \file KokkosBlas3_trmm.hpp

#include "KokkosKernels_Macros.hpp"
#include "KokkosBlas3_trmm_spec.hpp"
#include "KokkosKernels_helpers.hpp"
#include "KokkosKernels_Error.hpp"
#include <sstream>
#include <type_traits>

namespace KokkosBlas {

/// \brief Solve triangular linear system with multiple RHSs:
///        B = alpha * op(A) * B if side == "L" or "l"
///        B = alpha * B * op(A) if side == "R" or "r"
///
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
/// \tparam BViewType Input(RHS)/Output(solution) M-by-N matrix, as a 2-D
/// Kokkos::View
///
/// \param side  [in] "L" or "l" indicates matrix A is on the left of B
///                   "R" or "r" indicates matrix A is on the right of B
/// \param uplo  [in] "U" or "u" indicates matrix A is an upper triangular
/// matrix
///                   "L" or "l" indicates matrix A is a lower triangular matrix
/// \param trans [in] Specifies what op does to A:
//                    "N" or "n" for non-transpose,
//                    "T" or "t" for transpose,
//                    "C" or "c" for conjugate transpose.
/// \param diag  [in] "U" or "u" indicates the diagonal of A is assumed to be
/// unit
//                    "N" or "n" indicates the diagonal of A is assumed to be
//                    non-unit
/// \param alpha [in] Input coefficient used for
//                    multiplication with either A or B
/// \param A [in]     Input matrix, as a 2-D Kokkos::View
///                   If side == "L" or "l", matrix A is a M-by-M triangular
///                   matrix; otherwise, matrix A is a N-by-N triangular matrix
/// \param B [in,out] Input/Output matrix, as a 2-D Kokkos::View
///                   On entry, M-by-N matrix
///                   On exit, overwritten with the solution
template <class AViewType, class BViewType>
void trmm(const char side[], const char uplo[], const char trans[],
          const char diag[], typename BViewType::const_value_type& alpha,
          const AViewType& A, const BViewType& B) {
  static_assert(Kokkos::is_view<AViewType>::value,
                "AViewType must be a Kokkos::View.");
  static_assert(Kokkos::is_view<BViewType>::value,
                "BViewType must be a Kokkos::View.");
  static_assert(static_cast<int>(AViewType::rank) == 2,
                "AViewType must have rank 2.");
  static_assert(static_cast<int>(BViewType::rank) == 2,
                "BViewType must have rank 2.");

  // Check validity of indicator argument
  bool valid_side = (side[0] == 'L') || (side[0] == 'l') || (side[0] == 'R') ||
                    (side[0] == 'r');
  bool valid_uplo = (uplo[0] == 'U') || (uplo[0] == 'u') || (uplo[0] == 'L') ||
                    (uplo[0] == 'l');
  bool valid_trans = (trans[0] == 'N') || (trans[0] == 'n') ||
                     (trans[0] == 'T') || (trans[0] == 't') ||
                     (trans[0] == 'C') || (trans[0] == 'c');
  bool valid_diag = (diag[0] == 'U') || (diag[0] == 'u') || (diag[0] == 'N') ||
                    (diag[0] == 'n');
  if (!valid_side) {
    std::ostringstream os;
    os << "KokkosBlas::trmm: side = '" << side[0] << "'. "
       << "Valid values include 'L' or 'l' (A is on the left of X), "
          "'R' or 'r' (A is on the right of X).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (!valid_uplo) {
    std::ostringstream os;
    os << "KokkosBlas::trmm: uplo = '" << uplo[0] << "'. "
       << "Valid values include 'U' or 'u' (A is upper triangular), "
          "'L' or 'l' (A is lower triangular).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (!valid_trans) {
    std::ostringstream os;
    os << "KokkosBlas::trmm: trans = '" << trans[0] << "'. "
       << "Valid values include 'N' or 'n' (No transpose), 'T' or 't' "
          "(Transpose), "
          "and 'C' or 'c' (Conjugate transpose).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (!valid_diag) {
    std::ostringstream os;
    os << "KokkosBlas::trmm: diag = '" << diag[0] << "'. "
       << "Valid values include 'U' or 'u' (the diagonal of A is assumed to be "
          "unit), "
          "'N' or 'n' (the diagonal of A is assumed to be non-unit).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  bool is_A_lower_triangle = (side[0] == 'L' || side[0] == 'l');

  int64_t A_m = A.extent(0);
  int64_t A_n = A.extent(1);
  int64_t B_m = B.extent(0);
  int64_t B_n = B.extent(1);

  // Return if degenerated matrices are provided
  if (A_m == 0 || A_n == 0 || B_m == 0 || B_n == 0) return;

  // Ensure that the dimensions of A match and that we can legally perform A*B
  // or B*A
  if (A_m != A_n || (is_A_lower_triangle ? B_m : B_n) != A_n) {
    std::ostringstream os;
    os << "KokkosBlas::trmm: Dimensions of A and B do not match: "
       << "side: " << side[0] << " A: " << A.extent(0) << " x " << A.extent(1)
       << " B: " << B.extent(0) << " x " << B.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Create A matrix view type alias
  using AViewInternalType =
      Kokkos::View<typename AViewType::const_value_type**,
                   typename AViewType::array_layout,
                   typename AViewType::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  // Crease B matrix view type alias
  using BViewInternalType =
      Kokkos::View<typename BViewType::non_const_value_type**,
                   typename BViewType::array_layout,
                   typename BViewType::device_type,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  KokkosBlas::Impl::TRMM<AViewInternalType, BViewInternalType>::trmm(
      side, uplo, trans, diag, alpha, A, B);
}

}  // namespace KokkosBlas

#endif  // KOKKOS_BLAS3_TRMM_HPP_
