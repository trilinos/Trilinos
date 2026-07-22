// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_IDENTITY_HPP
#define KOKKOSBATCHED_IDENTITY_HPP

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"

#include "KokkosBatched_Copy_Decl.hpp"

namespace KokkosBatched {

/// \brief Batched Identity Operator:

class Identity {
 public:
  KOKKOS_INLINE_FUNCTION
  Identity() {}

  KOKKOS_INLINE_FUNCTION
  ~Identity() {}

  template <typename ArgTrans, typename ArgMode, int sameXY, typename MemberType, typename XViewType,
            typename YViewType>
  KOKKOS_INLINE_FUNCTION void apply(const MemberType &member, const XViewType &X, const YViewType &Y) const {
    if (sameXY == 0) {
      if (std::is_same<ArgMode, KokkosBatched::Mode::Serial>::value) {
        SerialCopy<Trans::NoTranspose>::invoke(X, Y);
      } else if (std::is_same<ArgMode, KokkosBatched::Mode::Team>::value) {
        TeamCopy<MemberType>::invoke(member, X, Y);
      } else if (std::is_same<ArgMode, KokkosBatched::Mode::TeamVector>::value) {
        TeamVectorCopy<MemberType>::invoke(member, X, Y);
      }
    }
  }
  template <typename ArgTrans, int sameXY, typename XViewType, typename YViewType>
  KOKKOS_INLINE_FUNCTION void apply(const XViewType &X, const YViewType &Y) const {
    if (sameXY == 0) {
      SerialCopy<Trans::NoTranspose>::invoke(X, Y);
    }
  }
};

}  // namespace KokkosBatched

#endif