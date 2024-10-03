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
#ifndef __KOKKOSBATCHED_HOUSEHOLDER_DECL_HPP__
#define __KOKKOSBATCHED_HOUSEHOLDER_DECL_HPP__

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
template <typename MemberType, typename ArgSide>
struct TeamVectorHouseholder {
  template <typename MemberType, typename aViewType, typename tauViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const aViewType &a, const tauViewType &tau);
};

}  // namespace KokkosBatched

#include "KokkosBatched_Householder_Serial_Impl.hpp"
#include "KokkosBatched_Householder_TeamVector_Impl.hpp"

#endif
