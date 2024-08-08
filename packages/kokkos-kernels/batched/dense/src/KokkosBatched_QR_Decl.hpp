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
#ifndef __KOKKOSBATCHED_QR_DECL_HPP__
#define __KOKKOSBATCHED_QR_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial QR
///

template <typename ArgAlgo>
struct SerialQR {
  template <typename AViewType, typename tViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const tViewType &t, const wViewType &w);
};

///
/// Team QR
///

template <typename MemberType, typename ArgAlgo>
struct TeamQR {
  template <typename AViewType, typename tViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType & /*member*/, const AViewType & /*A*/,
                                           const tViewType & /*t*/, const wViewType & /*w*/) {
    /// not implemented
    return -1;
  }
};

///
/// TeamVector QR
///

template <typename MemberType, typename ArgAlgo>
struct TeamVectorQR {
  template <typename AViewType, typename tViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const wViewType &w);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgMode, typename ArgAlgo>
struct QR {
  template <typename AViewType, typename tViewType, typename wViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                                const wViewType &w) {
    int r_val = 0;
    if (std::is_same<ArgMode, Mode::Serial>::value) {
      r_val = SerialQR<ArgAlgo>::invoke(A, t, w);
    } else if (std::is_same<ArgMode, Mode::Team>::value) {
      r_val = TeamQR<MemberType, ArgAlgo>::invoke(member, A, t, w);
    } else if (std::is_same<ArgMode, Mode::TeamVector>::value) {
      r_val = TeamVectorQR<MemberType, ArgAlgo>::invoke(member, A, t, w);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#include "KokkosBatched_QR_Serial_Impl.hpp"
#include "KokkosBatched_QR_TeamVector_Impl.hpp"

#endif
