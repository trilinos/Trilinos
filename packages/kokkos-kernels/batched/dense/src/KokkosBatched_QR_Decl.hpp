// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_QR_DECL_HPP
#define KOKKOSBATCHED_QR_DECL_HPP

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
    if constexpr (std::is_same_v<ArgMode, Mode::Serial>) {
      r_val = SerialQR<ArgAlgo>::invoke(A, t, w);
    } else if constexpr (std::is_same_v<ArgMode, Mode::Team>) {
      r_val = TeamQR<MemberType, ArgAlgo>::invoke(member, A, t, w);
    } else if constexpr (std::is_same_v<ArgMode, Mode::TeamVector>) {
      r_val = TeamVectorQR<MemberType, ArgAlgo>::invoke(member, A, t, w);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#include "KokkosBatched_QR_Serial_Impl.hpp"
#include "KokkosBatched_QR_TeamVector_Impl.hpp"

#endif
