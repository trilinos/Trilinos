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
#ifndef KOKKOSLAPACK_TRTRI_HPP_
#define KOKKOSLAPACK_TRTRI_HPP_

/// \file KokkosLapack_trtri.hpp

#include "KokkosKernels_Macros.hpp"
#include "KokkosLapack_trtri_spec.hpp"
#include "KokkosKernels_helpers.hpp"
#include <sstream>
#include <type_traits>
#include "KokkosKernels_Error.hpp"

namespace KokkosLapack {

/// \brief Find the inverse of the triangular matrix, A
///
///        A = inv(A)
///
/// \tparam AViewType Input matrix, as a 2-D Kokkos::View
///
/// \param uplo  [in] "U" or "u" indicates matrix A is an upper triangular
/// matrix
///                   "L" or "l" indicates matrix A is a lower triangular matrix
/// \param diag  [in] "U" or "u" indicates the diagonal of A is assumed to be
/// unit
//                    "N" or "n" indicates the diagonal of A is assumed to be
//                    non-unit
/// \param A [in,out] Input matrix, as a 2-D Kokkos::View
///                   On entry, A
///                   On successful exit, inv(A)
/// \return           0 upon success,
//                    i if the i-th diagonal elemet of A is zero, A is singular,
//                    and the inversion could not be completed.
// source: https://software.intel.com/en-us/mkl-developer-reference-c-trtri
template <class AViewType>
int trtri(const char uplo[], const char diag[], const AViewType& A) {
  static_assert(Kokkos::is_view<AViewType>::value, "AViewType must be a Kokkos::View.");
  static_assert(static_cast<int>(AViewType::rank) == 2, "AViewType must have rank 2.");

  // Check validity of indicator argument
  bool valid_uplo = (uplo[0] == 'U') || (uplo[0] == 'u') || (uplo[0] == 'L') || (uplo[0] == 'l');
  bool valid_diag = (diag[0] == 'U') || (diag[0] == 'u') || (diag[0] == 'N') || (diag[0] == 'n');

  if (!valid_uplo) {
    std::ostringstream os;
    os << "KokkosLapack::trtri: uplo = '" << uplo[0] << "'. "
       << "Valid values include 'U' or 'u' (A is upper triangular), "
          "'L' or 'l' (A is lower triangular).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }
  if (!valid_diag) {
    std::ostringstream os;
    os << "KokkosLapack::trtri: diag = '" << diag[0] << "'. "
       << "Valid values include 'U' or 'u' (the diagonal of A is assumed to be "
          "unit), "
          "'N' or 'n' (the diagonal of A is assumed to be non-unit).";
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  int64_t A_m = A.extent(0);
  int64_t A_n = A.extent(1);

  // Return if degenerated matrices are provided
  if (A_m == 0 || A_n == 0)
    return 0;  // This is success as the inverse of a matrix with no elements is
               // itself.

  // Ensure that the dimensions of A match and that we can legally perform A*B
  // or B*A
  if (A_m != A_n) {
    std::ostringstream os;
    os << "KokkosLapack::trtri: Dimensions of A do not match,"
       << " A: " << A.extent(0) << " x " << A.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Create A matrix view type alias
  using AViewInternalType = Kokkos::View<typename AViewType::non_const_value_type**, typename AViewType::array_layout,
                                         typename AViewType::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  // This is the return value type and should always reside on host
  using RViewInternalType =
      Kokkos::View<int, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  int result;
  RViewInternalType R = RViewInternalType(&result);

  KokkosLapack::Impl::TRTRI<RViewInternalType, AViewInternalType>::trtri(R, uplo, diag, A);

  return result;
}

}  // namespace KokkosLapack

#endif  // KOKKOSLAPACK_TRTRI_HPP_
