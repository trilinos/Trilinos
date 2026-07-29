// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMM_DECL_HPP
#define KOKKOSBATCHED_GEMM_DECL_HPP

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

/// \brief Serial Batched Gemm:
///
/// performs one of the matrix-matrix operations
///   C := alpha*op( A )*op( B ) + beta*C,
/// where op( X ) is one of
///   op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
///   alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,
///   op( B ) a k by n matrix and C an m by n matrix.
///
/// \tparam ArgTransA: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose) or A**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgTransB: Type indicating whether the B (Trans::NoTranspose), or B**T (Trans::Transpose) or B**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Gemm::Blocked) or unblocked
/// (KokkosBatched::Algo::Gemm::Unblocked) algorithm to be used
template <typename ArgTransA, typename ArgTransB, typename ArgAlgo>
struct SerialGemm {
  static_assert(KokkosBlas::is_trans_v<ArgTransA>, "KokkosBatched::SerialGemm: ArgTransA must be a KokkosBlas::Trans.");
  static_assert(KokkosBlas::is_trans_v<ArgTransB>, "KokkosBatched::SerialGemm: ArgTransB must be a KokkosBlas::Trans.");
#if defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL) && defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_BATCHED) && \
    defined(KOKKOSBATCHED_IMPL_ENABLE_INTEL_MKL_COMPACT_BATCHED)
  static_assert(
      std::is_same_v<ArgAlgo, Algo::Gemm::Unblocked> || std::is_same_v<ArgAlgo, Algo::Gemm::Blocked> ||
          std::is_same_v<ArgAlgo, Algo::Gemm::CompactMKL>,
      "KokkosBatched::SerialGemm: Use Algo::Gemm::Unblocked or Algo::Gemm::Blocked or Algo::Gemm::CompactMKL");
#else
  static_assert(std::is_same_v<ArgAlgo, Algo::Gemm::Unblocked> || std::is_same_v<ArgAlgo, Algo::Gemm::Blocked>,
                "KokkosBatched::SerialGemm: Use Algo::Gemm::Unblocked or Algo::Gemm::Blocked");
#endif

  /// \tparam ScalarType: Scalar type of alpha and beta
  /// \tparam AViewType: Input type for the matrix A, needs to be a 0D-2D view
  /// \tparam BViewType: Input type for the matrix B, needs to be a 0D-2D view
  /// \tparam CViewType: Input/Output type for the matrix C, needs to be a 0D-2D view
  ///
  /// \param[in] alpha: Scalar alpha
  /// \param[in] A: A is a dimension ( lda, ka ) matrix, where ka is k when ArgTransA = Trans::NoTranspose, and is m
  /// otherwise.
  /// \param[in] B: B is a dimension ( ldb, kb ) matrix, where kb is n when ArgTransB = Trans::NoTranspose, and is k
  /// otherwise.
  /// \param[in] beta: Scalar beta
  /// \param[in,out] C: C is a dimension ( ldc, n ) matrix. Before entry, the leading m by n part of the array C
  /// must contain the matrix C, except when beta is zero, in which case C need not be set on entry. On exit, the array
  /// C is overwritten by the m by n matrix ( alpha*op( A )*op( B ) + beta*C )
  ///
  /// No nested parallel_for is used inside of the function.
  ///
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType &A, const BViewType &B,
                                           const ScalarType beta, const CViewType &C);
};

/// \brief Team Batched Gemm:
///
/// performs one of the matrix-matrix operations
///   C := alpha*op( A )*op( B ) + beta*C,
/// where op( X ) is one of
///   op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
///   alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,
///   op( B ) a k by n matrix and C an m by n matrix.
///
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgTransA: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose) or A**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgTransB: Type indicating whether the B (Trans::NoTranspose), or B**T (Trans::Transpose) or B**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Gemm::Blocked) or unblocked
/// (KokkosBatched::Algo::Gemm::Unblocked) algorithm to be used
template <typename MemberType, typename ArgTransA, typename ArgTransB, typename ArgAlgo>
struct TeamGemm {
  static_assert(KokkosBlas::is_trans_v<ArgTransA>, "KokkosBatched::TeamGemm: ArgTransA must be a KokkosBlas::Trans.");
  static_assert(KokkosBlas::is_trans_v<ArgTransB>, "KokkosBatched::TeamGemm: ArgTransB must be a KokkosBlas::Trans.");
  static_assert(std::is_same_v<ArgAlgo, Algo::Gemm::Unblocked> || std::is_same_v<ArgAlgo, Algo::Gemm::Blocked>,
                "KokkosBatched::TeamGemm: Use Algo::Gemm::Unblocked or Algo::Gemm::Blocked");

  /// \tparam ScalarType: Scalar type of alpha and beta
  /// \tparam AViewType: Input type for the matrix A, needs to be a 0D-2D view
  /// \tparam BViewType: Input type for the matrix B, needs to be a 0D-2D view
  /// \tparam CViewType: Input/Output type for the matrix C, needs to be a 0D-2D view
  ///
  /// \param alpha [in]: Scalar alpha
  /// \param A [in]: A is a dimension ( lda, ka ) matrix, where ka is k when ArgTransA = Trans::NoTranspose, and is m
  /// otherwise.
  /// \param B [in]: B is a dimension ( ldb, kb ) matrix, where kb is n when ArgTransB = Trans::NoTranspose, and is k
  /// otherwise.
  /// \param beta [in]: Scalar beta
  /// \param C [inout]: C is a dimension ( ldc, n ) matrix. Before entry, the leading m by n part of the array C
  /// must contain the matrix C, except when beta is zero, in which case C need not be set on entry. On exit, the array
  /// C is overwritten by the m by n matrix ( alpha*op( A )*op( B ) + beta*C )
  ///
  /// Team thread parallelization is used inside of the function.
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C);
};

/// \brief TeamVector Batched Gemm:
///
/// performs one of the matrix-matrix operations
///   C := alpha*op( A )*op( B ) + beta*C,
/// where op( X ) is one of
///   op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
///   alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,
///   op( B ) a k by n matrix and C an m by n matrix.
/// \tparam MemberType: Member type of the Kokkos team policy
/// \tparam ArgTransA: Type indicating whether the A (Trans::NoTranspose), or A**T (Trans::Transpose) or A**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgTransB: Type indicating whether the B (Trans::NoTranspose), or B**T (Trans::Transpose) or B**H
/// (Trans::ConjTranspose) is used.
/// \tparam ArgAlgo: Type indicating the blocked (KokkosBatched::Algo::Gemm::Blocked) or unblocked
/// (KokkosBatched::Algo::Gemm::Unblocked) algorithm to be used
template <typename MemberType, typename ArgTransA, typename ArgTransB, typename ArgAlgo>
struct TeamVectorGemm {
  static_assert(KokkosBlas::is_trans_v<ArgTransA>,
                "KokkosBatched::TeamVectorGemm: ArgTransA must be a KokkosBlas::Trans.");
  static_assert(KokkosBlas::is_trans_v<ArgTransB>,
                "KokkosBatched::TeamVectorGemm: ArgTransB must be a KokkosBlas::Trans.");
  static_assert(std::is_same_v<ArgAlgo, Algo::Gemm::Unblocked> || std::is_same_v<ArgAlgo, Algo::Gemm::Blocked>,
                "KokkosBatched::Gemm: Use Algo::Gemm::Unblocked or Algo::Gemm::Blocked");

  /// \tparam ScalarType: Scalar type of alpha and beta
  /// \tparam AViewType: Input type for the matrix A, needs to be a 0D-2D view
  /// \tparam BViewType: Input type for the matrix B, needs to be a 0D-2D view
  /// \tparam CViewType: Input/Output type for the matrix C, needs to be a 0D-2D view
  ///
  /// \param alpha [in]: Scalar alpha
  /// \param A [in]: A is a dimension ( lda, ka ) matrix, where ka is k when ArgTransA = Trans::NoTranspose, and is m
  /// otherwise.
  /// \param B [in]: B is a dimension ( ldb, kb ) matrix, where kb is n when ArgTransB = Trans::NoTranspose, and is k
  /// otherwise.
  /// \param beta [in]: Scalar beta
  /// \param C [inout]: C is a dimension ( ldc, n ) matrix. Before entry, the leading m by n part of the array C
  /// must contain the matrix C, except when beta is zero, in which case C need not be set on entry. On exit, the array
  /// C is overwritten by the m by n matrix ( alpha*op( A )*op( B ) + beta*C )
  ///
  /// Team vector parallelization is used inside of the function.
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgTransA, typename ArgTransB, typename ArgMode, typename ArgAlgo>
struct Gemm {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                                const BViewType &B, const ScalarType beta, const CViewType &C) {
    int r_val = 0;
    if constexpr (std::is_same_v<ArgMode, Mode::Serial>) {
      r_val = SerialGemm<ArgTransA, ArgTransB, ArgAlgo>::invoke(alpha, A, B, beta, C);
    } else if constexpr (std::is_same_v<ArgMode, Mode::Team>) {
      r_val = TeamGemm<MemberType, ArgTransA, ArgTransB, ArgAlgo>::invoke(member, alpha, A, B, beta, C);
    } else if constexpr (std::is_same_v<ArgMode, Mode::TeamVector>) {
      r_val = TeamVectorGemm<MemberType, ArgTransA, ArgTransB, ArgAlgo>::invoke(member, alpha, A, B, beta, C);
    }
    return r_val;
  }
};
}  // namespace KokkosBatched

#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_Gemm_Team_Impl.hpp"
#include "KokkosBatched_Gemm_TeamVector_Impl.hpp"

#endif  // KOKKOSBATCHED_GEMM_DECL_HPP
