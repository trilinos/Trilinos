// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMM_TEAM_IMPL_HPP
#define KOKKOSBATCHED_GEMM_TEAM_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemm_Team_Internal.hpp"

namespace KokkosBatched {

///
/// Team Impl
/// =========

///
/// Implemented:
/// NT/NT, T/NT, NT/T, T/T
///
/// Not yet implemented (ConjTranspose)
/// CT/NT, NT/CT, CT/CT
///

///
/// NT/NT
///

template <typename MemberType>
struct TeamGemm<MemberType, Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_1 = Impl::get_extent_int(A, 1);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Unblocked>::invoke(member, C_extent_0, C_extent_1, A_extent_1, alpha, A.data(),
                                                           A_stride_0, A_stride_1, B.data(), B_stride_0, B_stride_1,
                                                           beta, C.data(), C_stride_0, C_stride_1);
  }
};

template <typename MemberType>
struct TeamGemm<MemberType, Trans::NoTranspose, Trans::NoTranspose, Algo::Gemm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_1 = Impl::get_extent_int(A, 1);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Blocked>::invoke(member, C_extent_0, C_extent_1, A_extent_1, alpha, A.data(),
                                                         A_stride_0, A_stride_1, B.data(), B_stride_0, B_stride_1, beta,
                                                         C.data(), C_stride_0, C_stride_1);
  }
};

///
/// T/NT
///

template <typename MemberType>
struct TeamGemm<MemberType, Trans::Transpose, Trans::NoTranspose, Algo::Gemm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Unblocked>::invoke(member, C_extent_0, C_extent_1, A_extent_0, alpha, A.data(),
                                                           A_stride_1, A_stride_0, B.data(), B_stride_0, B_stride_1,
                                                           beta, C.data(), C_stride_0, C_stride_1);
  }
};

template <typename MemberType>
struct TeamGemm<MemberType, Trans::Transpose, Trans::NoTranspose, Algo::Gemm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Blocked>::invoke(member, C_extent_0, C_extent_1, A_extent_0, alpha, A.data(),
                                                         A_stride_1, A_stride_0, B.data(), B_stride_0, B_stride_1, beta,
                                                         C.data(), C_stride_0, C_stride_1);
  }
};

///
/// NT/T
///

template <typename MemberType>
struct TeamGemm<MemberType, Trans::NoTranspose, Trans::Transpose, Algo::Gemm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_1 = Impl::get_extent_int(A, 1);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Unblocked>::invoke(member, C_extent_0, C_extent_1, A_extent_1, alpha, A.data(),
                                                           A_stride_0, A_stride_1, B.data(), B_stride_1, B_stride_0,
                                                           beta, C.data(), C_stride_0, C_stride_1);
  }
};

template <typename MemberType>
struct TeamGemm<MemberType, Trans::NoTranspose, Trans::Transpose, Algo::Gemm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_1 = Impl::get_extent_int(A, 1);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Blocked>::invoke(member, C_extent_0, C_extent_1, A_extent_1, alpha, A.data(),
                                                         A_stride_0, A_stride_1, B.data(), B_stride_1, B_stride_0, beta,
                                                         C.data(), C_stride_0, C_stride_1);
  }
};

///
/// T/T
///

template <typename MemberType>
struct TeamGemm<MemberType, Trans::Transpose, Trans::Transpose, Algo::Gemm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Unblocked>::invoke(member, C_extent_0, C_extent_1, A_extent_0, alpha, A.data(),
                                                           A_stride_1, A_stride_0, B.data(), B_stride_1, B_stride_0,
                                                           beta, C.data(), C_stride_0, C_stride_1);
  }
};

template <typename MemberType>
struct TeamGemm<MemberType, Trans::Transpose, Trans::Transpose, Algo::Gemm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return TeamGemmInternal<Algo::Gemm::Blocked>::invoke(member, C_extent_0, C_extent_1, A_extent_0, alpha, A.data(),
                                                         A_stride_1, A_stride_0, B.data(), B_stride_1, B_stride_0, beta,
                                                         C.data(), C_stride_0, C_stride_1);
  }
};

}  // namespace KokkosBatched

#endif
