// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_GEMM_TEAM_IMPL_HPP
#define KOKKOSBATCHED_GEMM_TEAM_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)
/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBlas_util.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Gemm_Common_Impl.hpp"
#include "KokkosBatched_Gemm_Team_Internal.hpp"

namespace KokkosBatched {

///
/// Team Impl
/// =========

///
/// NT/NT
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::NoTranspose, Trans::NoTranspose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_1 = Impl::get_extent_int(A, 1);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_1;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::NoTranspose, Trans::NoTranspose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A B
    // C (m x n), A(m x k), B(k x n)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C_extent_0, C_extent_1, A_extent_1, alpha, A.data(),
        A_stride_0, A_stride_1, B.data(), B_stride_0, B_stride_1, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// T/NT
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::Transpose, Trans::NoTranspose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_0;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::Transpose, Trans::NoTranspose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A^T B
    // C (m x n), A(k x m), B(k x n)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C_extent_0, C_extent_1, A_extent_0, alpha, A.data(),
        A_stride_1, A_stride_0, B.data(), B_stride_0, B_stride_1, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// C/NT
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::ConjTranspose, Trans::NoTranspose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_0;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::ConjTranspose, Trans::NoTranspose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A^H B
    // C (m x n), A(k x m), B(k x n)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpID(), C_extent_0, C_extent_1, A_extent_0, alpha,
        A.data(), A_stride_1, A_stride_0, B.data(), B_stride_0, B_stride_1, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// NT/T
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::NoTranspose, Trans::Transpose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_1 = Impl::get_extent_int(A, 1);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_1;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::NoTranspose, Trans::Transpose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A B^T
    // C (m x n), A(m x k), B(n x k)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C_extent_0, C_extent_1, A_extent_1, alpha, A.data(),
        A_stride_0, A_stride_1, B.data(), B_stride_1, B_stride_0, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// T/T
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::Transpose, Trans::Transpose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_0;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::Transpose, Trans::Transpose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A^T B^T
    // C (m x n), A(k x m), B(n x k)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpID(), C_extent_0, C_extent_1, A_extent_0, alpha, A.data(),
        A_stride_1, A_stride_0, B.data(), B_stride_1, B_stride_0, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// C/T
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::ConjTranspose, Trans::Transpose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_0;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::ConjTranspose, Trans::Transpose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A^H B^T
    // C (m x n), A(k x m), B(n x k)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpID(), C_extent_0, C_extent_1, A_extent_0, alpha,
        A.data(), A_stride_1, A_stride_0, B.data(), B_stride_1, B_stride_0, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// NT/C
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::NoTranspose, Trans::ConjTranspose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_1 = Impl::get_extent_int(A, 1);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_1;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::NoTranspose, Trans::ConjTranspose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A B^H
    // C (m x n), A(m x k), B(n x k)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpConj(), C_extent_0, C_extent_1, A_extent_1, alpha,
        A.data(), A_stride_0, A_stride_1, B.data(), B_stride_1, B_stride_0, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// T/C
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::Transpose, Trans::ConjTranspose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_0;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::Transpose, Trans::ConjTranspose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A^T B^H
    // C (m x n), A(k x m), B(n x k)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpID(), KokkosBlas::Impl::OpConj(), C_extent_0, C_extent_1, A_extent_0, alpha,
        A.data(), A_stride_1, A_stride_0, B.data(), B_stride_1, B_stride_0, beta, C.data(), C_stride_0, C_stride_1);
  }
};

///
/// C/C
///

template <typename MemberType, typename ArgAlgo>
struct TeamGemm<MemberType, Trans::ConjTranspose, Trans::ConjTranspose, ArgAlgo> {
  template <typename ScalarType, typename AViewType, typename BViewType, typename CViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B, const ScalarType beta, const CViewType &C) {
    const int A_extent_0 = Impl::get_extent_int(A, 0);
    const int C_extent_0 = Impl::get_extent_int(C, 0), C_extent_1 = Impl::get_extent_int(C, 1);

    const std::size_t A_stride_0 = Impl::get_stride(A, 0), A_stride_1 = Impl::get_stride(A, 1);
    const std::size_t B_stride_0 = Impl::get_stride(B, 0), B_stride_1 = Impl::get_stride(B, 1);
    const std::size_t C_stride_0 = Impl::get_stride(C, 0), C_stride_1 = Impl::get_stride(C, 1);

    // Quick return if possible
    const int m = C_extent_0, n = C_extent_1, k = A_extent_0;
    if (m == 0 || n == 0 || ((alpha == ScalarType(0) || k == 0) && beta == ScalarType(1))) return 0;

    auto info = Impl::checkGemmInput<Trans::ConjTranspose, Trans::ConjTranspose>(A, B, C);
    if (info) return info;

    // C = beta C + alpha A^H B^H
    // C (m x n), A(k x m), B(n x k)
    return Impl::TeamGemmInternal<ArgAlgo>::invoke(
        member, KokkosBlas::Impl::OpConj(), KokkosBlas::Impl::OpConj(), C_extent_0, C_extent_1, A_extent_0, alpha,
        A.data(), A_stride_1, A_stride_0, B.data(), B_stride_1, B_stride_0, beta, C.data(), C_stride_0, C_stride_1);
  }
};

}  // namespace KokkosBatched

#endif
