// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SYRK_HPP_
#define KOKKOSBATCHED_SYRK_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Syrk:
/// Performs one of the symmetric rank k operation
///   C := alpha*A*A**T + beta*C
///   C := alpha*A*A**H + beta*C
///   C := alpha*A**T*A + beta*C
///   C := alpha*A**H*A + beta*C
///    where alpha and beta are real scalars, C is an n by n symmetric or hermitian matrix and
///    A is an n by k matrix in the first case and a k by n matrix in the second case.
///
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose), or A**H
/// (Trans::ConjTranspose), or AH (Trans::ConjNoTranspose) is used.
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgUplo, typename ArgTrans>
struct SerialSyrk {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syrk: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::syrk: Use Trans::NoTranspose/Trans::Transpose for {s,d,c,z}syrk or "
                "Trans::ConjNoTranspose/Trans::ConjTranspose for {c,z}herk");

  template <typename ScalarType, typename AViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const ScalarType beta,
                                           const CViewType &C);
};

/// \brief Team Batched Syrk:
/// Performs one of the symmetric rank k operation
///   C := alpha*A*A**T + beta*C
///   C := alpha*A*A**H + beta*C
///   C := alpha*A**T*A + beta*C
///   C := alpha*A**H*A + beta*C
///    where alpha and beta are real scalars, C is an n by n hermitian matrix and
///    A is an n by k matrix in the first case and a k by n matrix in the second case.
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose), or A**H
/// (Trans::ConjTranspose), or AH (Trans::ConjNoTranspose) is used.
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamSyrk {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syrk: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::syrk: Use Trans::NoTranspose/Trans::Transpose for {s,d,c,z}syrk or "
                "Trans::ConjNoTranspose/Trans::ConjTranspose for {c,z}herk");
  template <typename ScalarType, typename AViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const ScalarType beta, const CViewType &C);
};

/// \brief TeamVector Batched Syrk:
/// Performs one of the symmetric rank k operation
///   C := alpha*A*A**T + beta*C
///   C := alpha*A*A**H + beta*C
///   C := alpha*A**T*A + beta*C
///   C := alpha*A**H*A + beta*C
///    where alpha and beta are real scalars, C is an n by n hermitian matrix and
///    A is an n by k matrix in the first case and a k by n matrix in the second case.
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is
/// modified
/// \tparam ArgTrans: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose), or A**H
/// (Trans::ConjTranspose), or AH (Trans::ConjNoTranspose) is used.
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamVectorSyrk {
  static_assert(
      std::is_same_v<ArgUplo, Uplo::Upper> || std::is_same_v<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syrk: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::syrk: Use Trans::NoTranspose/Trans::Transpose for {s,d,c,z}syrk or "
                "Trans::ConjNoTranspose/Trans::ConjTranspose for {c,z}herk");
  template <typename ScalarType, typename AViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const ScalarType beta, const CViewType &C);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Syrk_Impl.hpp"

#endif  // KOKKOSBATCHED_SYRK_HPP_
