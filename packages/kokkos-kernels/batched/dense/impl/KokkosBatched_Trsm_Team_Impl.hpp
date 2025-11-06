// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_TRSM_TEAM_IMPL_HPP
#define KOKKOSBATCHED_TRSM_TEAM_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trsm_Team_Internal.hpp"

namespace KokkosBatched {

///
/// Team Impl
/// =========

///
/// L/L/NT
///
/// B := inv(tril(A)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B_extent_1, alpha, A.data(), A.stride(0),
                                                                    A.stride(1), B.data(), B.stride(0), B_stride_1);
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B_extent_1, alpha, A.data(), A.stride(0), A.stride(1),
                                                                  B.data(), B.stride(0), B_stride_1);
  }
};

///
/// R/U/NT
///
/// B := (alpha*B) inv(triu(A))
/// A(n x n), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B_extent_1,
                                                                    B.extent(0), alpha, A.data(), A.stride(1),
                                                                    A.stride(0), B.data(), B_stride_1, B.stride(0));
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B_extent_1,
                                                                  B.extent(0), alpha, A.data(), A.stride(1),
                                                                  A.stride(0), B.data(), B_stride_1, B.stride(0));
  }
};

///
/// R/L/NT
///
/// B := (alpha*B) inv(tril(A))
/// A(n x n), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B_extent_1,
                                                                    B.extent(0), alpha, A.data(), A.stride(1),
                                                                    A.stride(0), B.data(), B_stride_1, B.stride(0));
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Lower, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B_extent_1,
                                                                  B.extent(0), alpha, A.data(), A.stride(1),
                                                                  A.stride(0), B.data(), B_stride_1, B.stride(0));
  }
};

///
/// R/U/T
///
/// B := (alpha*B) inv(triu(A))
/// A(n x n), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B_extent_1,
                                                                    B.extent(0), alpha, A.data(), A.stride(0),
                                                                    A.stride(1), B.data(), B_stride_1, B.stride(0));
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Right, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B_extent_1,
                                                                  B.extent(0), alpha, A.data(), A.stride(0),
                                                                  A.stride(1), B.data(), B_stride_1, B.stride(0));
  }
};

///
/// L/U/NT
///
/// B := inv(triu(A)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B_extent_1, alpha, A.data(), A.stride(0),
                                                                    A.stride(1), B.data(), B.stride(0), B_stride_1);
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::NoTranspose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B_extent_1, alpha, A.data(), A.stride(0), A.stride(1),
                                                                  B.data(), B.stride(0), B_stride_1);
  }
};

///
/// L/L/T
///
/// B := inv(tril(AT)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B_extent_1, alpha, A.data(), A.stride(1),
                                                                    A.stride(0), B.data(), B.stride(0), B_stride_1);
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Lower, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftUpper<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B_extent_1, alpha, A.data(), A.stride(1), A.stride(0),
                                                                  B.data(), B.stride(0), B_stride_1);
  }
};

///
/// L/U/T
///
/// B := inv(triu(AT)) (alpha*B)
/// A(m x m), B(m x n)

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Unblocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftLower<Algo::Trsm::Unblocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                    B_extent_1, alpha, A.data(), A.stride(1),
                                                                    A.stride(0), B.data(), B.stride(0), B_stride_1);
  }
};

template <typename MemberType, typename ArgDiag>
struct TeamTrsm<MemberType, Side::Left, Uplo::Upper, Trans::Transpose, ArgDiag, Algo::Trsm::Blocked> {
  template <typename ScalarType, typename AViewType, typename BViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const ScalarType alpha, const AViewType &A,
                                           const BViewType &B) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamTrsmInternalLeftLower<Algo::Trsm::Blocked>::invoke(member, ArgDiag::use_unit_diag, B.extent(0),
                                                                  B_extent_1, alpha, A.data(), A.stride(1), A.stride(0),
                                                                  B.data(), B.stride(0), B_stride_1);
  }
};

}  // namespace KokkosBatched

#endif
