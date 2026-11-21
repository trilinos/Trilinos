// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SPR_HPP_
#define KOKKOSBATCHED_SPR_HPP_

#include <concepts>
#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Batched Serial Spr:
/// Performs the symmetric rank 1 operation
///   A := alpha*x*x**T + A or A := alpha*x*x**H + A
///    where alpha is a scalar, x is an n element vector, and A is a n by n symmetric or Hermitian matrix, supplied in
///    packed form.
///
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the transpose (Trans::Transpose) or conjugate transpose
/// (Trans::ConjTranspose) of x is used
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \tparam APViewType: Input/output type for the matrix A, needs to be a 1D view (packed storage)
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] x: x is a length n vector, a rank 1 view
/// \param[inout] ap: ap is a n by n matrix, supplied in packed form, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgUplo, typename ArgTrans>
struct SerialSpr {
  static_assert(
      std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
      "KokkosBatched::spr: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::spr: Use Trans::Transpose for {s,d,c,z}spr or Trans::ConjTranspose for {c,z}her");
  template <typename ScalarType, typename XViewType, typename APViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const XViewType &x, const APViewType &ap);
};

/// \brief Batched Team Spr:
/// Performs the symmetric rank 1 operation
///   A := alpha*x*x**T + A or A := alpha*x*x**H + A
///    where alpha is a scalar, x is an n element vector, and A is a n by n symmetric or Hermitian matrix, supplied in
///    packed form.
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the transpose (Trans::Transpose) or conjugate transpose
/// (Trans::ConjTranspose) of x is used
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \tparam APViewType: Input/output type for the matrix A, needs to be a 1D view (packed storage)
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] x: x is a length n vector, a rank 1 view
/// \param[inout] ap: ap is a n by n matrix, supplied in packed form, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamSpr {
  static_assert(
      std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
      "KokkosBatched::spr: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::spr: Use Trans::Transpose for {s,d,c,z}spr or Trans::ConjTranspose for {c,z}her");
  template <typename ScalarType, typename XViewType, typename APViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const XViewType &x,
                                           const APViewType &ap);
};

/// \brief Batched TeamVector Spr:
/// Performs the symmetric rank 1 operation
///   A := alpha*x*x**T + A or A := alpha*x*x**H + A
///    where alpha is a scalar, x is an n element vector, and A is a n by n symmetric or Hermitian matrix, supplied in
///    packed form.
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the transpose (Trans::Transpose) or conjugate transpose
/// (Trans::ConjTranspose) of x is used
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam XViewType: Input type for the vector x, needs to be a 1D view
/// \tparam APViewType: Input/output type for the matrix A, needs to be a 1D view (packed storage)
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] x: x is a length n vector, a rank 1 view
/// \param[inout] ap: ap is a n by n matrix, supplied in packed form, a rank 1 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamVectorSpr {
  static_assert(
      std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
      "KokkosBatched::spr: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::spr: Use Trans::Transpose for {s,d,c,z}spr or Trans::ConjTranspose for {c,z}her");
  template <typename ScalarType, typename XViewType, typename APViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const XViewType &x,
                                           const APViewType &ap);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Spr_Impl.hpp"

#endif  // KOKKOSBATCHED_SPR_HPP_
