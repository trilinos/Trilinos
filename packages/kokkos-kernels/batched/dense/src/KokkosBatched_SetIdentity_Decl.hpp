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
#ifndef __KOKKOSBATCHED_SET_IDENTITY_DECL_HPP__
#define __KOKKOSBATCHED_SET_IDENTITY_DECL_HPP__

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
