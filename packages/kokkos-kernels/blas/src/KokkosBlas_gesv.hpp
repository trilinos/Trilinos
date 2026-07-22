// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosBlas_gesv.hpp
/// \brief Local dense linear solve
///
/// This file provides KokkosBlas::gesv. This function performs a
/// local (no MPI) dense linear solve on a system of linear equations
/// A * X = B where A is a general N-by-N matrix and X and B are N-by-NRHS
/// matrices.

#ifndef KOKKOSBLAS_GESV_HPP_
#define KOKKOSBLAS_GESV_HPP_

#include "KokkosLapack_gesv.hpp"

namespace KokkosBlas {

/// \brief Solve the dense linear equation system A*X = B.
///
/// \tparam AMatrix Input matrix/Output LU, as a 2-D Kokkos::View.
/// \tparam BXMV Input (right-hand side)/Output (solution) (multi)vector, as a
/// 1-D or 2-D Kokkos::View. \tparam IPIVV Output pivot indices, as a 1-D
/// Kokkos::View
///
/// \param A [in,out] On entry, the N-by-N matrix to be solved. On exit, the
/// factors L and U from
///   the factorization A = P*L*U; the unit diagonal elements of L are not
///   stored.
/// \param B [in,out] On entry, the right hand side (multi)vector B. On exit,
/// the solution (multi)vector X. \param IPIV [out] On exit, the pivot indices
/// (for partial pivoting). If the View extents are zero and
///   its data pointer is NULL, pivoting is not used.
///
template <class AMatrix, class BXMV, class IPIVV>
[[deprecated]] void gesv(const AMatrix& A, const BXMV& B, const IPIVV& IPIV) {
  KokkosLapack::gesv(A, B, IPIV);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS_GESV_HPP_
