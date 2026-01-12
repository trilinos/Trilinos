// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_HOUSEHOLDER_DECL_HPP
#define KOKKOSBATCHED_HOUSEHOLDER_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Householder
///

// level 1 operation
template <typename ArgSide>
struct SerialHouseholder {
  template <typename aViewType, typename tauViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const aViewType &a, const tauViewType &tau);
};

///
/// TeamVector Householder
///

// level 1 operation
template <typename ArgSide>
struct TeamVectorHouseholder {
  template <typename MemberType, typename aViewType, typename tauViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const aViewType &a, const tauViewType &tau);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Householder_Serial_Impl.hpp"
#include "KokkosBatched_Householder_TeamVector_Impl.hpp"

#endif
