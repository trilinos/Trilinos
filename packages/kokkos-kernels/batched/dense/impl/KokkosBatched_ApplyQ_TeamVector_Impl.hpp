//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef __KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_APPLY_Q_TEAMVECTOR_IMPL_HPP__

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
    return TeamVectorApplyQ_LeftForwardInternal::invoke(member, B.extent(0), B.extent(1), A.extent(1), A.data(),
                                                        A.stride_0(), A.stride_1(), t.data(), t.stride_0(), B.data(),
                                                        B.stride_0(), B.stride_1(), w.data());
  }
};

template <typename MemberType>
struct TeamVectorApplyQ<MemberType, Side::Left, Trans::Transpose, Algo::ApplyQ::Unblocked> {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w) {
    return TeamVectorApplyQ_LeftBackwardInternal::invoke(member, B.extent(0), B.extent(1), A.extent(1), A.data(),
                                                         A.stride_0(), A.stride_1(), t.data(), t.stride_0(), B.data(),
                                                         B.stride_0(), B.stride_1(), w.data());
  }
};

template <typename MemberType>
struct TeamVectorApplyQ<MemberType, Side::Right, Trans::NoTranspose, Algo::ApplyQ::Unblocked> {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w) {
    return TeamVectorApplyQ_RightForwardInternal::invoke(member, B.extent(0), B.extent(1), A.extent(1), A.data(),
                                                         A.stride_0(), A.stride_1(), t.data(), t.stride_0(), B.data(),
                                                         B.stride_0(), B.stride_1(), w.data());
  }
};

}  // namespace KokkosBatched

#endif
