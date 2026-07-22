// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_APPLY_Q_DECL_HPP
#define KOKKOSBATCHED_APPLY_Q_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial ApplyQ
///

template <typename ArgSide, typename ArgTrans, typename ArgAlgo>
struct SerialApplyQ {
  static_assert(is_side_v<ArgSide>, "KokkosBatched::SerialApplyQ: ArgSide must be a KokkosBatched::Side.");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::SerialApplyQ: ArgTrans must be a KokkosBlas::Trans.");
  static_assert(KokkosBlas::is_level2_v<ArgAlgo>,
                "KokkosBatched::SerialApplyQ: ArgAlgo must be a KokkosBlas::Algo::Level2.");

  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const AViewType &A, const tViewType &t, const BViewType &B,
                                           const wViewType &w);
};

///
/// Team ApplyQ
///

template <typename MemberType, typename ArgSide, typename ArgTrans, typename ArgAlgo>
struct TeamApplyQ {
  static_assert(is_side_v<ArgSide>, "KokkosBatched::SerialApplyQ: ArgSide must be a KokkosBatched::Side.");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::SerialApplyQ: ArgTrans must be a KokkosBlas::Trans.");
  static_assert(KokkosBlas::is_level2_v<ArgAlgo>,
                "KokkosBatched::SerialApplyQ: ArgAlgo must be a KokkosBlas::Algo::Level2.");

  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w);
};

///
/// TeamVector ApplyQ
///

template <typename MemberType, typename ArgSide, typename ArgTrans, typename ArgAlgo>
struct TeamVectorApplyQ {
  static_assert(is_side_v<ArgSide>, "KokkosBatched::SerialApplyQ: ArgSide must be a KokkosBatched::Side.");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::SerialApplyQ: ArgTrans must be a KokkosBlas::Trans.");
  static_assert(KokkosBlas::is_level2_v<ArgAlgo>,
                "KokkosBatched::SerialApplyQ: ArgAlgo must be a KokkosBlas::Algo::Level2.");

  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                           const BViewType &B, const wViewType &w);
};

///
/// Selective Interface
///
template <typename MemberType, typename ArgSide, typename ArgTrans, typename ArgMode, typename ArgAlgo>
struct ApplyQ {
  static_assert(is_side_v<ArgSide>, "KokkosBatched::ApplyQ: ArgSide must be a KokkosBatched::Side.");
  static_assert(KokkosBlas::is_trans_v<ArgTrans>, "KokkosBatched::ApplyQ: ArgTrans must be a KokkosBlas::Trans.");
  static_assert(KokkosBlas::is_mode_v<ArgMode>, "KokkosBatched::ApplyQ: ArgMode must be a KokkosBlas::Mode.");
  static_assert(KokkosBlas::is_level2_v<ArgAlgo>, "KokkosBatched::ApplyQ: ArgAlgo must be a KokkosBlas::Algo::Level2.");

  template <typename AViewType, typename tViewType, typename BViewType, typename wViewType>
  KOKKOS_FORCEINLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const tViewType &t,
                                                const BViewType &B, const wViewType &w) {
    int r_val = 0;
    if constexpr (std::is_same_v<ArgMode, Mode::Serial>) {
      r_val = SerialApplyQ<ArgSide, ArgTrans, ArgAlgo>::invoke(A, t, B, w);
    } else if constexpr (std::is_same_v<ArgMode, Mode::Team>) {
      r_val = TeamApplyQ<MemberType, ArgSide, ArgTrans, ArgAlgo>::invoke(member, A, t, B, w);
    } else if constexpr (std::is_same_v<ArgMode, Mode::TeamVector>) {
      r_val = TeamVectorApplyQ<MemberType, ArgSide, ArgTrans, ArgAlgo>::invoke(member, A, t, B, w);
    }
    return r_val;
  }
};

}  // namespace KokkosBatched

#include "KokkosBatched_ApplyQ_Serial_Impl.hpp"
#include "KokkosBatched_ApplyQ_TeamVector_Impl.hpp"

#endif
