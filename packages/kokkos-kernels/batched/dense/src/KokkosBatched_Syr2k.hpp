// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SYR2K_HPP_
#define KOKKOSBATCHED_SYR2K_HPP_

#include <concepts>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

/// \brief Serial Batched Syr2k:
/// Performs one of the symmetric rank 2k operations
///   C := alpha*A*B^T + alpha*B*A^T + beta*C
///   C := alpha*A*B^H + conj(alpha)*B*A^H + beta*C
///   C := alpha*A^T*B + alpha*B^T*A + beta*C
///   C := alpha*A^H*B + conj(alpha)*B^H*A + beta*C
///    where alpha and beta are real scalars, C is an n by n symmetric or hermitian matrix and
///    A and B are n by k matrices in the first case and k by n matrices in the second case.
///
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of C is
/// modified
/// \tparam ArgTrans: Type indicating the operation to be performed as follows:
///   ArgTrans::NoTranspose: C := alpha*A*B^T + alpha*B*A^T + beta*C
///   ArgTrans::ConjNoTranspose: C := alpha*A*B^H + conj(alpha)*B*A^H + beta*C
///   ArgTrans::Transpose: C := alpha*A^T*B + alpha*B^T*A + beta*C
///   ArgTrans::ConjTranspose: C := alpha*A^H*B + conj(alpha)*B^H*A + beta*C
///
/// \tparam AlphaType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam BViewType: Input type for the matrix B, needs to be a 2D view
/// \tparam BetaType: Input type for the scalar beta
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] B: B is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] beta: beta is a scalar
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
/// No nested parallel_for is used inside of the function.
///
template <typename ArgUplo, typename ArgTrans>
struct SerialSyr2k {
  static_assert(
      std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syr2k: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::syr2k: Use Trans::NoTranspose/Trans::Transpose for {s,d,c,z}syr2k or "
                "Trans::ConjNoTranspose/Trans::ConjTranspose for {c,z}her2k");

  template <typename AlphaType, typename AViewType, typename BViewType, typename BetaType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AlphaType alpha, const AViewType &A, const BViewType &B,
                                           const BetaType beta, const CViewType &C);
};

/// \brief Team Batched Syr2k:
/// Performs one of the symmetric rank 2k operations
///   C := alpha*A*B^T + alpha*B*A^T + beta*C
///   C := alpha*A*B^H + conj(alpha)*B*A^H + beta*C
///   C := alpha*A^T*B + alpha*B^T*A + beta*C
///   C := alpha*A^H*B + alpha*B^H*A + beta*C
///    where alpha and beta are real scalars, C is an n by n symmetric or hermitian matrix and
///    A and B are n by k matrices in the first case and k by n matrices in the second case.
///
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of C is
/// modified
/// \tparam ArgTrans: Type indicating the operation to be performed as follows:
///   ArgTrans::NoTranspose: C := alpha*A*B^T + alpha*B*A^T + beta*C
///   ArgTrans::ConjNoTranspose: C := alpha*A*B^H + conj(alpha)*B*A^H + beta*C
///   ArgTrans::Transpose: C := alpha*A^T*B + alpha*B^T*A + beta*C
///   ArgTrans::ConjTranspose: C := alpha*A^H*B + conj(alpha)*B^H*A + beta*C
///
/// \tparam AlphaType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam BViewType: Input type for the matrix B, needs to be a 2D view
/// \tparam BetaType: Input type for the scalar beta
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] member: Kokkos team member handle
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] B: B is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] beta: beta is a scalar
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
/// Team thread parallelization is used inside of the function.
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamSyr2k {
  static_assert(
      std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syr2k: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::syr2k: Use Trans::NoTranspose/Trans::Transpose for {s,d,c,z}syr2k or "
                "Trans::ConjNoTranspose/Trans::ConjTranspose for {c,z}her2k");
  template <typename AlphaType, typename AViewType, typename BViewType, typename BetaType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AlphaType alpha, const AViewType &A,
                                           const BViewType &B, const BetaType beta, const CViewType &C);
};

/// \brief TeamVector Batched Syr2k:
/// Performs one of the symmetric rank 2k operations
///   C := alpha*A*B^T + alpha*B*A^T + beta*C
///   C := alpha*A*B^H + conj(alpha)*B*A^H + beta*C
///   C := alpha*A^T*B + alpha*B^T*A + beta*C
///   C := alpha*A^H*B + conj(alpha)*B^H*A + beta*C
///    where alpha and beta are real scalars, C is an n by n symmetric or hermitian matrix and
///    A and B are n by k matrices in the first case and k by n matrices in the second case.
///
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of C is
/// modified
/// \tparam ArgTrans: Type indicating the operation to be performed as follows:
///   ArgTrans::NoTranspose: C := alpha*A*B^T + alpha*B*A^T + beta*C
///   ArgTrans::ConjNoTranspose: C := alpha*A*B^H + conj(alpha)*B*A^H + beta*C
///   ArgTrans::Transpose: C := alpha*A^T*B + alpha*B^T*A + beta*C
///   ArgTrans::ConjTranspose: C := alpha*A^H*B + conj(alpha)*B^H*A + beta*C
///
/// \tparam AlphaType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam BViewType: Input type for the matrix B, needs to be a 2D view
/// \tparam BetaType: Input type for the scalar beta
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] member: Kokkos team member handle
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] B: B is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] beta: beta is a scalar
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
/// Team thread parallelization is used inside of the function.
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamVectorSyr2k {
  static_assert(
      std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
      "KokkosBatched::syr2k: Use Uplo::Upper for upper triangular matrix or Uplo::Lower for lower triangular matrix");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>,
                "KokkosBatched::syr2k: Use Trans::NoTranspose/Trans::Transpose for {s,d,c,z}syr2k or "
                "Trans::ConjNoTranspose/Trans::ConjTranspose for {c,z}her2k");
  template <typename AlphaType, typename AViewType, typename BViewType, typename BetaType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AlphaType alpha, const AViewType &A,
                                           const BViewType &B, const BetaType beta, const CViewType &C);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Syr2k_Impl.hpp"

#endif  // KOKKOSBATCHED_SYR2K_HPP_
