// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SYMV_HPP
#define KOKKOSBATCHED_SYMV_HPP

#include <concepts>
#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

/// \brief Serial Batched Symv:
///
/// performs one of the matrix-vector operations
///   y := alpha*A*x + beta*y,
///   alpha and beta are scalars, x and y are n element vectors, and A is an n by n symmetric or hermitian matrix.
///
/// \tparam ArgUplo: Type indicating whether the A is upper (Uplo::Upper) or lower (Uplo::Lower) triangular.
/// \tparam ArgTrans: Type indicating whether A is symmetric (Trans::Transpose) or hermitian (Trans::ConjTranspose)
/// matrix.
template <typename ArgUplo, typename ArgTrans>
struct SerialSymv {
  static_assert(std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
                "KokkosBatched::symv: ArgUplo must be a KokkosBatched::Uplo.");
  static_assert(std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::symv: Use Trans::Transpose for {s,d,c,z}symv or Trans::ConjTranspose for {c,z}hemv");

  /// \tparam ScalarType: Scalar type of alpha and beta
  /// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
  /// \tparam XViewType: Input type for the vector x, needs to be a 1D view
  /// \tparam YViewType: Input/Output type for the vector y, needs to be a 1D view
  ///
  /// \param[in] alpha: Scalar alpha
  /// \param[in] A: A is a dimension ( lda, n ) matrix. Before entry with uplo = Uplo::Upper, the leading n by n upper
  /// triangular part of the array A must contain the upper triangular part of the symmetric (hermitian) matrix and the
  /// strictly lower triangular part of A is not referenced. Before entry with uplo = Uplo::Lower, the leading n by n
  /// lower triangular part of the array A must contain the lower triangular part of the symmetric (hermitian) matrix
  /// and the strictly upper triangular part of A is not referenced.
  /// \param[in] x: x is a dimension ( n ) vector
  /// \param[in] beta: Scalar beta
  /// \param[in,out] y: y is a dimension ( n ) vector. Before entry, y must contain the vector y. On exit, y is
  /// overwritten by the result ( alpha*A*x + beta*y )
  ///
  /// No nested parallel_for is used inside of the function.
  ///
  template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const XViewType &x,
                                           const ScalarType beta, const YViewType &y);
};

/// \brief Team Batched Symv:
///
/// performs one of the matrix-vector operations
///   y := alpha*A*x + beta*y,
///   alpha and beta are scalars, x and y are n element vectors, and A is an n by n symmetric or hermitian matrix.
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgUplo: Type indicating whether the A is upper (Uplo::Upper) or lower (Uplo::Lower) triangular.
/// \tparam ArgTrans: Type indicating whether A is symmetric (Trans::Transpose) or hermitian (Trans::ConjTranspose)
/// matrix.
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamSymv {
  static_assert(std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
                "KokkosBatched::symv: ArgUplo must be a KokkosBatched::Uplo.");
  static_assert(std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::symv: Use Trans::Transpose for {s,d,c,z}symv or Trans::ConjTranspose for {c,z}hemv");

  /// \tparam ScalarType: Scalar type of alpha and beta
  /// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
  /// \tparam XViewType: Input type for the vector x, needs to be a 1D view
  /// \tparam YViewType: Input/Output type for the vector y, needs to be a 1D view
  ///
  /// \param[in] alpha: Scalar alpha
  /// \param[in] A: A is a dimension ( lda, n ) matrix. Before entry with uplo = Uplo::Upper, the leading n by n upper
  /// triangular part of the array A must contain the upper triangular part of the symmetric (hermitian) matrix and the
  /// strictly lower triangular part of A is not referenced. Before entry with uplo = Uplo::Lower, the leading n by n
  /// lower triangular part of the array A must contain the lower triangular part of the symmetric (hermitian) matrix
  /// and the strictly upper triangular part of A is not referenced.
  /// \param[in] x: x is a dimension ( n ) vector
  /// \param[in] beta: Scalar beta
  /// \param[in,out] y: y is a dimension ( n ) vector. Before entry, y must contain the vector y. On exit, y is
  /// overwritten by the result ( alpha*A*x + beta*y )
  ///
  /// Team thread parallelization is used inside of the function.
  ///
  template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const XViewType &x, const ScalarType beta, const YViewType &y);
};

/// \brief TeamVector Batched Symv:
///
/// performs one of the matrix-vector operations
///   y := alpha*A*x + beta*y,
///   alpha and beta are scalars, x and y are n element vectors, and A is an n by n symmetric or hermitian matrix.
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgUplo: Type indicating whether the A is upper (Uplo::Upper) or lower (Uplo::Lower) triangular.
/// \tparam ArgTrans: Type indicating whether A is symmetric (Trans::Transpose) or hermitian (Trans::ConjTranspose)
/// matrix.
template <typename MemberType, typename ArgUplo, typename ArgTrans>
struct TeamVectorSymv {
  static_assert(std::same_as<ArgUplo, Uplo::Upper> || std::same_as<ArgUplo, Uplo::Lower>,
                "KokkosBatched::symv: ArgUplo must be a KokkosBatched::Uplo.");
  static_assert(std::same_as<ArgTrans, Trans::Transpose> || std::same_as<ArgTrans, Trans::ConjTranspose>,
                "KokkosBatched::symv: Use Trans::Transpose for {s,d,c,z}symv or Trans::ConjTranspose for {c,z}hemv");

  /// \tparam ScalarType: Scalar type of alpha and beta
  /// \tparam AViewType: Input type for the matrix A, needs to be a 2D view
  /// \tparam XViewType: Input type for the vector x, needs to be a 1D view
  /// \tparam YViewType: Input/Output type for the vector y, needs to be a 1D view
  ///
  /// \param[in] alpha: Scalar alpha
  /// \param[in] A: A is a dimension ( lda, n ) matrix. Before entry with uplo = Uplo::Upper, the leading n by n upper
  /// triangular part of the array A must contain the upper triangular part of the symmetric (hermitian) matrix and the
  /// strictly lower triangular part of A is not referenced. Before entry with uplo = Uplo::Lower, the leading n by n
  /// lower triangular part of the array A must contain the lower triangular part of the symmetric (hermitian) matrix
  /// and the strictly upper triangular part of A is not referenced.
  /// \param[in] x: x is a dimension ( n ) vector
  /// \param[in] beta: Scalar beta
  /// \param[in,out] y: y is a dimension ( n ) vector. Before entry, y must contain the vector y. On exit, y is
  /// overwritten by the result ( alpha*A*x + beta*y )
  ///
  /// Team vector parallelization is used inside of the function.
  ///
  template <typename ScalarType, typename AViewType, typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const XViewType &x, const ScalarType beta, const YViewType &y);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Symv_Impl.hpp"

#endif  // KOKKOSBATCHED_SYMV_HPP
