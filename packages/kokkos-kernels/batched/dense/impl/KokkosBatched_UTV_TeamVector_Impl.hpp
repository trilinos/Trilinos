// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_UTV_TEAMVECTOR_IMPL_HPP
#define KOKKOSBATCHED_UTV_TEAMVECTOR_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_UTV_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
struct TeamVectorUTV<MemberType, Algo::UTV::Unblocked> {
  template <typename AViewType, typename pViewType, typename UViewType, typename VViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const pViewType &p,
                                           const UViewType &U, const VViewType &V, const wViewType &w,
                                           int &matrix_rank) {
    return TeamVectorUTV_Internal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1),
                                          p.data(), p.stride(0), U.data(), U.stride(0), U.stride(1), V.data(),
                                          V.stride(0), V.stride(1), w.data(), matrix_rank);
  }
};

}  // namespace KokkosBatched

#endif
