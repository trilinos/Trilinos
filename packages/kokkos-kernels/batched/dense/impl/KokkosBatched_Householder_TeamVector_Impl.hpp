// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_HOUSEHOLDER_TEAMVECTOR_IMPL_HPP
#define KOKKOSBATCHED_HOUSEHOLDER_TEAMVECTOR_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <>
template <typename MemberType, typename aViewType, typename tauViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorHouseholder<Side::Left>::invoke(const MemberType &member, const aViewType &a,
                                                                     const tauViewType &tau) {
  return TeamVectorLeftHouseholderInternal::invoke(member, a.extent(0) - 1, a.data(), a.data() + a.stride(0),
                                                   a.stride(0), tau.data());
}

}  // namespace KokkosBatched

#endif
