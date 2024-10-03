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
#ifndef __KOKKOSBATCHED_APPLY_Q_DECL_HPP__
#define __KOKKOSBATCHED_APPLY_Q_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial ApplyQ
///

template <typename ArgSide, typename ArgTrans, typename ArgAlgo>
struct SerialApplyQ {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const tViewType &t, const BViewType &B,
                                           const wViewType &w);
};

///
/// Team ApplyQ
///

template <typename MemberType, typename ArgSide, typename ArgTrans, typename ArgAlgo>
struct TeamApplyQ {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w);
};

///
/// TeamVector ApplyQ
///

template <typename MemberType, typename ArgSide, typename ArgTrans, typename ArgAlgo>
struct TeamVectorApplyQ {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgSide, typename ArgTrans, typename ArgMode, typename ArgAlgo>
struct ApplyQ {
  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                                const BViewType &B, const wViewType &w) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialApplyQ<ArgSide, ArgTrans, ArgAlgo>::invoke(A, t, B, w);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamApplyQ<MemberType, ArgSide, ArgTrans, ArgAlgo>::invoke(member, A, t, B, w);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamVectorApplyQ<MemberType, ArgSide, ArgTrans, ArgAlgo>::invoke(member, A, t, B, w);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#include "KokkosBatched_ApplyQ_Serial_Impl.hpp"
#include "KokkosBatched_ApplyQ_TeamVector_Impl.hpp"

#endif
