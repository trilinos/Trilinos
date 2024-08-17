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
#ifndef __KOKKOSBATCHED_HOUSEHOLDER_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_HOUSEHOLDER_TEAMVECTOR_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
template <typename aViewType, typename tauViewType>
KOKKOS_INLINE_FUNCTION int TeamVectorHouseholder<Side::Left>::invoke(const MemberType &member, const aViewType &a,
                                                                     const tauViewType &tau) {
  return TeamVectorLeftHouseholderInternal::invoke(member, a.extent(0) - 1, a.data(), a.data() + a.stride(0),
                                                   a.stride(0), tau.data());
}

}  // namespace KokkosBatched

#endif
