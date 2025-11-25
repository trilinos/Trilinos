// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP
#define KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Householder
///

// level 1 operation
template <typename ArgSide, typename ArgTrans = Trans::NoTranspose>
struct SerialApplyHouseholder {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const uViewType &u2, const tauViewType &tau, const AViewType &A,
                                           const wViewType &w);
};

// level 1 operation
template <typename MemberType, typename ArgSide>
struct TeamVectorApplyHouseholder {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const uViewType &u2, const tauViewType &tau,
                                           const AViewType &A, const wViewType &w);
};

}  // namespace KokkosBatched

#include "KokkosBatched_ApplyHouseholder_Serial_Impl.hpp"
#include "KokkosBatched_ApplyHouseholder_TeamVector_Impl.hpp"

#endif
