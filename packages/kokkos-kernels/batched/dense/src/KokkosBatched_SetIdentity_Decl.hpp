// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_SET_IDENTITY_DECL_HPP
#define KOKKOSBATCHED_SET_IDENTITY_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"

namespace KokkosBatched {

///
/// Serial SetIdentity
///

struct SerialSetIdentity {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A);
};

///
/// Team Set
///

template <typename MemberType>
struct TeamSetIdentity {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgMode>
struct SetIdentity {
  template <typename AViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialSetIdentity::invoke(A);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamSetIdentity<MemberType>::invoke(member, A);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#endif
