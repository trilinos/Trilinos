// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_TRMV_HPP_
#define KOKKOSBATCHED_TRMV_HPP_

#include <concepts>
#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Trmv:
/// Performs the matrix-vector operation
///   x := alpha*op( A )*x
/// where op( A ) is either A, A**T, or A**H depending on the value of ArgTrans.
/// alpha is a scalar, x is an n element vector, and A is a n by n matrix.
///
/// \tparam ArgUplo: Type indicating whether the matrix is an upper (Uplo::Upper) or lower (Uplo::Lower) triangular
/// matrix modified
/// \tparam ArgTrans: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose) or A**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgDiag: Type indicating whether the matrix is unit triangular (Diag::Unit) or non-unit triangular
/// (Diag::NonUnit)
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \tparam AViewType: Input/output type for the matrix A, needs to be a 2D view
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by n matrix, a rank 2 view
/// \param[in,out] x: x is a length n vector, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgUplo, typename ArgTrans, typename ArgDiag>
struct SerialTrmv {
  static_assert(
      is_uplo_v<ArgUplo>,
      "KokkosBatched::trmv: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::trmv: ArgTrans must be a KokkosBlas::Trans");
  static_assert(is_diag_v<ArgDiag>, "KokkosBatched::trmv: ArgDiag must be either Diag::Unit or Diag::NonUnit");
  template <typename ScalarType, typename AViewType, typename XViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const XViewType &x);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Trmv_Impl.hpp"

#endif  // KOKKOSBATCHED_TRMV_HPP_
