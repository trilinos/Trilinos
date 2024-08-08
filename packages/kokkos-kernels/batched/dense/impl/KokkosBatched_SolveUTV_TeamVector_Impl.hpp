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
#ifndef __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_IMPL_HPP__
#define __KOKKOSBATCHED_SOLVE_UTV_TEAMVECTOR_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_SolveUTV_TeamVector_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Impl
/// ===============
template <typename MemberType>
struct TeamVectorSolveUTV<MemberType, Algo::UTV::Unblocked> {
  template <typename UViewType, typename TViewType, typename VViewType, typename pViewType, typename XViewType,
            typename BViewType, typename wViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int matrix_rank, const UViewType &U,
                                           const TViewType &T, const VViewType &V, const pViewType &p,
                                           const XViewType &X, const BViewType &B, const wViewType &w) {
    if (BViewType::rank == 1)
      TeamVectorSolveUTV_Internal::invoke(member, matrix_rank, T.extent(0), V.extent(0), U.data(), U.stride(0),
                                          U.stride(1), T.data(), T.stride(0), T.stride(1), V.data(), V.stride(0),
                                          V.stride(1), p.data(), p.stride(0), X.data(), X.stride(0), B.data(),
                                          B.stride(0), w.data());
    else
      TeamVectorSolveUTV_Internal::invoke(member, matrix_rank, T.extent(0), V.extent(0), B.extent(1), U.data(),
                                          U.stride(0), U.stride(1), T.data(), T.stride(0), T.stride(1), V.data(),
                                          V.stride(0), V.stride(1), p.data(), p.stride(0), X.data(), X.stride(0),
                                          X.stride(1), B.data(), B.stride(0), B.stride(1), w.data());
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
