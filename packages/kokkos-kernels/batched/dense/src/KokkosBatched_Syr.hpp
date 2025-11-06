// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SYR_HPP_
#define KOKKOSBATCHED_SYR_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Syr:
/// Performs the symmetric rank 1 operation
///   A := alpha*x*x**T + A or A := alpha*x*x**H + A
///    where alpha is a scalar, x is an n element vector, and A is a n by n symmetric or Hermitian matrix.
///
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the transpose (Trans::Transpose) or conjugate transpose
/// (Trans::ConjTranspose) of x is used
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \tparam AViewType: Input/output type for the matrix A, needs to be a 2D view
///
/// \param alpha [in]: alpha is a scalar
/// \param x [in]: x is a length n vector, a rank 1 view
/// \param A [inout]: A is a n by n matrix, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgUplo, typename ArgTrans>
struct SerialSyr {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syr: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::is_same_v<ArgTrans, Trans::Transpose> || std::is_same_v<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::syr: Use Trans::Transpose for {s,d,c,z}syr or Trans::ConjTranspose for {c,z}her");
  template <typename ScalarType, typename XViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const AViewType &a);
};
}  // namespace KokkosBatched

#include "KokkosBatched_Syr_Serial_Impl.hpp"

#endif  // KOKKOSBATCHED_SYR_HPP_
