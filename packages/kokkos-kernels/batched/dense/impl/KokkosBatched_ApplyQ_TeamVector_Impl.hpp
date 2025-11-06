// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_IMPL_HPP
#define KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_ApplyQ_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
struct TeamVectorApplyQ<MemberType, Side::Left, Trans::NoTranspose, Algo::ApplyQ::Unblocked> {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamVectorApplyQ_LeftForwardInternal::invoke(member, B.extent(0), B_extent_1, A.extent(1), A.data(),
                                                        A.stride(0), A.stride(1), t.data(), t.stride(0), B.data(),
                                                        B.stride(0), B_stride_1, w.data());
  }
};

template <typename MemberType>
struct TeamVectorApplyQ<MemberType, Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked> {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamVectorApplyQ_LeftBackwardInternal::invoke(member, B.extent(0), B_extent_1, A.extent(1), A.data(),
                                                         A.stride(0), A.stride(1), t.data(), t.stride(0), B.data(),
                                                         B.stride(0), B_stride_1, w.data());
  }
};

template <typename MemberType>
struct TeamVectorApplyQ<MemberType, Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked> {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w) {
    static_assert(AViewType::rank() == 2);
    constexpr size_t B_rank = BViewType::rank();
    static_assert(B_rank == 1 || B_rank == 2);

    // Quick return if possible
    if (B.size() == 0) return 0;

    size_t B_extent_1 = B_rank == 1 ? 1 : B.extent(1);
    size_t B_stride_1 = B_rank == 1 ? 1 : B.stride(1);

    return TeamVectorApplyQ_RightForwardInternal::invoke(member, B.extent(0), B_extent_1, A.extent(1), A.data(),
                                                         A.stride(0), A.stride(1), t.data(), t.stride(0), B.data(),
                                                         B.stride(0), B_stride_1, w.data());
  }
};

}  // namespace KokkosBatched

#endif
