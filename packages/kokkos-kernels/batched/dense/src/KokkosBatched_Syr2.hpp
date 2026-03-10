// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SYR2_HPP_
#define KOKKOSBATCHED_SYR2_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Syr2:
/// Performs the symmetric rank 2 operation
///   A := alpha*x*y**T + alpha*y*x**T + A or A := alpha*x*y**H + conjg(alpha)*y*x**H + A
///    where alpha is a scalar, x and y are n element vectors, and A is a n by n symmetric or Hermitian matrix.
///
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the transpose (Trans::Transpose) or conjugate transpose
/// (Trans::ConjTranspose) of x and y is used
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \tparam YViewType: Input type for the vector y, needs to be a 1D view
/// \tparam AViewType: Input/output type for the matrix A, needs to be a 2D view
///
/// \param[in] alpha : alpha is a scalar
/// \param[in] x : x is a length n vector, a rank 1 view
/// \param[in] y : y is a length n vector, a rank 1 view
/// \param[inout] A : A is a n by n matrix, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgUplo, typename ArgTrans>
struct SerialSyr2 {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syr2: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::is_same_v<ArgTrans, Trans::Transpose> || std::is_same_v<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::syr2: Use Trans::Transpose for {s,d,c,z}syr2 or Trans::ConjTranspose for {c,z}her2");
  template <typename ScalarType, typename XViewType, typename YViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const YViewType &y,
                                           const AViewType &a);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Syr2_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_SYR2_HPP_
