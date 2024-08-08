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
#ifndef __KOKKOSBATCHED_UTV_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_UTV_TEAMVECTOR_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_UTV_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============

template <typename MemberType>
struct TeamVectorUTV<MemberType, Algo::UTV::Unblocked> {
  template <typename AViewType, typename pViewType, typename UViewType, typename VViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const AViewType &A, const pViewType &p,
                                           const UViewType &U, const VViewType &V, const wViewType &w,
                                           int &matrix_rank) {
    return TeamVectorUTV_Internal::invoke(member, A.extent(0), A.extent(1), A.data(), A.stride(0), A.stride(1),
                                          p.data(), p.stride(0), U.data(), U.stride(0), U.stride(1), V.data(),
                                          V.stride(0), V.stride(1), w.data(), matrix_rank);
  }
};

}  // namespace KokkosBatched

#endif
