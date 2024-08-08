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
#ifndef KOKKOSBLAS_TRTRI_HPP_
#define KOKKOSBLAS_TRTRI_HPP_

/// \file KokkosBlas_trtri.hpp

#include "KokkosLapack_trtri.hpp"

namespace KokkosBlas {

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
[[deprecated]] int trtri(const char uplo[], const char diag[], const AViewType& A) {
  return KokkosLapack::trtri(uplo, diag, A);
}

}  // namespace KokkosBlas

#endif  // KOKKOS_BLASLAPACK_TRTRI_HPP_
