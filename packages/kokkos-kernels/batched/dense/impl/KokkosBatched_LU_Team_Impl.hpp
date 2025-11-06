// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_LU_TEAM_IMPL_HPP
#define KOKKOSBATCHED_LU_TEAM_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_LU_Team_Internal.hpp"

namespace KokkosBatched {

///
/// Team Impl
/// =========

///
/// LU no piv
///

template <typename MemberType>
struct TeamLU<MemberType, Algo::LU::Unblocked> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const AViewType &A,
      const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0) {
    return TeamLU_Internal<Algo::LU::Unblocked>::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride(0),
                                                        A.stride(1), tiny);
  }
};

template <typename MemberType>
struct TeamLU<MemberType, Algo::LU::Blocked> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const MemberType &member, const AViewType &A,
      const typename MagnitudeScalarType<typename AViewType::non_const_value_type>::type tiny = 0) {
    return TeamLU_Internal<Algo::LU::Blocked>::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride(0),
                                                      A.stride(1), tiny);
  }
};

}  // namespace KokkosBatched

#endif
