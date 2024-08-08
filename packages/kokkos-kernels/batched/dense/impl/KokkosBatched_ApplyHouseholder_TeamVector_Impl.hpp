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
#ifndef __KOKKOSBATCHED_APPLY_HOUSEHOLDER_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_APPLY_HOUSEHOLDER_TEAMVECTOR_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Householder_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// Serial Impl
/// ===========

template <typename MemberType>
struct TeamVectorApplyHouseholder<MemberType, Side::Left> {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const uViewType &u2, const tauViewType &tau,
                                           const AViewType &A, const wViewType &w) {
    return TeamVectorApplyLeftHouseholderInternal::invoke(member, A.extent(0) - 1, A.extent(1), tau.data(), u2.data(),
                                                          u2.stride(0), A.data(), A.stride(1), A.data() + A.stride(0),
                                                          A.stride(0), A.stride(1), w.data());
  }
};

template <typename MemberType>
struct TeamVectorApplyHouseholder<MemberType, Side::Right> {
  template <typename uViewType, typename tauViewType, typename AViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const uViewType &u2, const tauViewType &tau,
                                           const AViewType &A, const wViewType &w) {
    return TeamVectorApplyRightHouseholderInternal::invoke(member, A.extent(0), A.extent(1) - 1, tau.data(), u2.data(),
                                                           u2.stride(0), A.data(), A.stride(0), A.data() + A.stride(1),
                                                           A.stride(0), A.stride(1), w.data());
  }
};

}  // namespace KokkosBatched

#endif
