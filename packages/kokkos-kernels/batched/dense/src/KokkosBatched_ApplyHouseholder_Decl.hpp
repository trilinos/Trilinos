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
#ifndef __KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP__
#define __KOKKOSBATCHED_APPLY_HOUSEHOLDER_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Householder
///

// level 1 operation
template <typename ArgSide>
struct SerialApplyHouseholder {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const uViewType &u2, const tauViewType &tau,
                                           const AViewType const wViewType &w);
};

// level 1 operation
template <typename MemberType, typename ArgSide>
struct TeamVectorApplyHouseholder {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const uViewType &u2, const tauViewType &tau,
                                           const AViewType const wViewType &w);
};

}  // namespace KokkosBatched

#include "KokkosBatched_ApplyHouseholder_Serial_Impl.hpp"
#include "KokkosBatched_ApplyHouseholder_TeamVector_Impl.hpp"

#endif
