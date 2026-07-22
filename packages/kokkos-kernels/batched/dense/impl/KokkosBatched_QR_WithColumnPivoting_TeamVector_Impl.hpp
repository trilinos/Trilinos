// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_QR_WITH_COLUMNPIVOTING_TEAMVECTOR_IMPL_HPP
#define KOKKOSBATCHED_QR_WITH_COLUMNPIVOTING_TEAMVECTOR_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_QR_WithColumnPivoting_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
struct TeamVectorQR_WithColumnPivoting<MemberType, Algo::QR::Unblocked> {
  template <typename AViewType, typename tViewType, typename pViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const pViewType &p, const wViewType &w,
                                           /* */ int &matrix_rank) {
    return TeamVectorQR_WithColumnPivotingInternal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride(0),
                                                           A.stride(1), t.data(), t.stride(0), p.data(), p.stride(0),
                                                           w.data(), matrix_rank);
  }
};

}  // namespace KokkosBatched

#endif
