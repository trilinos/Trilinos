// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_SYMM_HPP_
#define KOKKOSBATCHED_SYMM_HPP_

#include <concepts>
#include <Kokkos_Core.hpp>
#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
/// \brief Serial Batched Symm:
/// Performs one of the matrix-matrix operations
///   C := alpha*A*B + beta*C
///   C := alpha*B*A + beta*C
///    where alpha and beta are real scalars, A is a symmetric or hermitian matrix and
///    B and C are m by n matrices
///
/// \tparam ArgSide: Type indicating whether the symmetric or hermitian matrix A appears on the left (Side::Left) or
/// right (Side::Right) in the operation
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is to
/// be referenced modified
/// \tparam Conj: Boolean indicating whether the matrix A is symmetric (false) or hermitian (true)
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam BViewType: Input type for the matrix B, needs to be a 2D view
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
template <typename ArgSide, typename ArgUplo, bool Conj>
struct SerialSymm {
  static_assert(is_side_v<ArgSide>, "KokkosBatched::symm: ArgSide must be a KokkosBatched::Side.");
  static_assert(is_uplo_v<ArgUplo>, "KokkosBatched::symm: ArgUplo must be a KokkosBatched::Uplo.");
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B,
                                           const ScalarType beta, const CViewType &C);
};

/// \brief Team Batched Symm:
/// Performs one of the matrix-matrix operations
///   C := alpha*A*B + beta*C
///   C := alpha*B*A + beta*C
///    where alpha and beta are real scalars, A is a symmetric or hermitian matrix and
///    B and C are m by n matrices
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgSide: Type indicating whether the symmetric or hermitian matrix A appears on the left (Side::Left) or
/// right (Side::Right) in the operation
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is to
/// be referenced modified
/// \tparam Conj: Boolean indicating whether the matrix A is symmetric (false) or hermitian (true)
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam BViewType: Input type for the matrix B, needs to be a 2D view
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] B: B is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] beta: beta is a scalar
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
/// Team thread parallelization is used inside of the function.
///
template <typename MemberType, typename ArgSide, typename ArgUplo, bool Conj>
struct TeamSymm {
  static_assert(is_side_v<ArgSide>, "KokkosBatched::symm: ArgSide must be a KokkosBatched::Side.");
  static_assert(is_uplo_v<ArgUplo>, "KokkosBatched::symm: ArgUplo must be a KokkosBatched::Uplo.");
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C);
};

/// \brief TeamVector Batched Symm:
/// Performs one of the matrix-matrix operations
///   C := alpha*A*B + beta*C
///   C := alpha*B*A + beta*C
///    where alpha and beta are real scalars, A is a symmetric or hermitian matrix and
///    B and C are m by n matrices
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgSide: Type indicating whether the symmetric or hermitian matrix A appears on the left (Side::Left) or
/// right (Side::Right) in the operation
/// \tparam ArgUplo: Type indicating whether the upper (Uplo::Upper) or lower (Uplo::Lower) triangular part of A is to
/// be referenced modified
/// \tparam Conj: Boolean indicating whether the matrix A is symmetric (false) or hermitian (true)
///
/// \tparam ScalarType: Input type for the scalar alpha
/// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
/// \tparam BViewType: Input type for the matrix B, needs to be a 2D view
/// \tparam CViewType: Input/output type for the matrix C, needs to be a 2D view
///
/// \param[in] alpha: alpha is a scalar
/// \param[in] A: A is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] B: B is a n by k matrix, a rank 2 view (or k by n if ArgTrans indicates transpose)
/// \param[in] beta: beta is a scalar
/// \param[in,out] C: C is a n by n matrix, a rank 2 view
///
/// Team vector parallelization is used inside of the function.
///
template <typename MemberType, typename ArgSide, typename ArgUplo, bool Conj>
struct TeamVectorSymm {
  static_assert(is_side_v<ArgSide>, "KokkosBatched::symm: ArgSide must be a KokkosBatched::Side.");
  static_assert(is_uplo_v<ArgUplo>, "KokkosBatched::symm: ArgUplo must be a KokkosBatched::Uplo.");
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Symm_Impl.hpp"

#endif  // KOKKOSBATCHED_SYMM_HPP_
