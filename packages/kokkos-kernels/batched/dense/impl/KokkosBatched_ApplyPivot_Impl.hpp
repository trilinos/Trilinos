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
#ifndef __KOKKOSBATCHED_APPLY_PIVOT_IMPL_HPP__
#define __KOKKOSBATCHED_APPLY_PIVOT_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_ApplyPivot_Internal.hpp"

namespace KokkosBatched {

///
/// TeamVector Internal Impl
/// ========================

///
/// Forward pivot apply
///

/// row swap
template <typename MemberType>
struct TeamVectorApplyPivot<MemberType, Side::Left, Direct::Forward> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int as0 = A.stride(0);
      TeamVectorApplyPivotVectorForwardInternal::invoke(member, piv, A.data(), as0);
    } else if (AViewType::rank == 2) {
      const int n = A.extent(1), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixForwardInternal::invoke(member, n, piv, A.data(), as0, as1);
    }
    return 0;
  }

  template <typename PivViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const PivViewType piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int plen = piv.extent(0), ps0 = piv.stride(0), as0 = A.stride(0);
      TeamVectorApplyPivotVectorForwardInternal::invoke(member, plen, piv.data(), ps0, A.data(), as0);
    } else if (AViewType::rank == 2) {
      // row permutation
      const int plen = piv.extent(0), ps0 = piv.stride(0), n = A.extent(1), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixForwardInternal::invoke(member, n, plen, piv.data(), ps0, A.data(), as0, as1);
    }
    return 0;
  }
};

/// column swap
template <typename MemberType>
struct TeamVectorApplyPivot<MemberType, Side::Right, Direct::Forward> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int as0 = A.stride(0);
      TeamVectorApplyPivotVectorForwardInternal::invoke(member, piv, A.data(), as0);
    } else if (AViewType::rank == 2) {
      const int m = A.extent(0), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixForwardInternal::invoke(member, m, piv, A.data(), as1, as0);
    }
    return 0;
  }

  template <typename PivViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const PivViewType &piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int plen = piv.extent(0), as0 = A.stride(0);
      TeamVectorApplyPivotVectorForwardInternal ::invoke(member, plen, piv.data(), A.data(), as0);
    } else if (AViewType::rank == 2) {
      // column permutation
      const int plen = piv.extent(0), ps = piv.stride(0), m = A.extent(0), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixForwardInternal ::invoke(member, m, plen, piv.data(), ps, A.data(), as1, as0);
    }
    return 0;
  }
};

///
/// Backward pivot apply
///

/// row swap
template <typename MemberType>
struct TeamVectorApplyPivot<MemberType, Side::Left, Direct::Backward> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int as0 = A.stride(0);
      TeamVectorApplyPivotVectorBackwardInternal::invoke(member, piv, A.data(), as0);
    } else if (AViewType::rank == 2) {
      const int n = A.extent(1), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixBackwardInternal::invoke(member, n, piv, A.data(), as0, as1);
    }
    return 0;
  }

  template <typename PivViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const PivViewType piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int plen = piv.extent(0), ps0 = piv.stride(0), as0 = A.stride(0);
      TeamVectorApplyPivotVectorBackwardInternal::invoke(member, plen, piv.data(), ps0, A.data(), as0);
    } else if (AViewType::rank == 2) {
      // row permutation
      const int plen = piv.extent(0), ps0 = piv.stride(0), n = A.extent(1), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixBackwardInternal::invoke(member, n, plen, piv.data(), ps0, A.data(), as0, as1);
    }
    return 0;
  }
};

/// column swap
template <typename MemberType>
struct TeamVectorApplyPivot<MemberType, Side::Right, Direct::Backward> {
  template <typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int as0 = A.stride(0);
      TeamVectorApplyPivotVectorBackwardInternal::invoke(member, piv, A.data(), as0);
    } else if (AViewType::rank == 2) {
      const int m = A.extent(0), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixBackwardInternal::invoke(member, m, piv, A.data(), as1, as0);
    }
    return 0;
  }

  template <typename PivViewType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const PivViewType &piv, const AViewType &A) {
    if (AViewType::rank == 1) {
      const int plen = piv.extent(0), as0 = A.stride(0);
      TeamVectorApplyPivotVectorBackwardInternal ::invoke(member, plen, piv.data(), A.data(), as0);
    } else if (AViewType::rank == 2) {
      // column permutation
      const int plen = piv.extent(0), ps = piv.stride(0), m = A.extent(0), as0 = A.stride(0), as1 = A.stride(1);
      TeamVectorApplyPivotMatrixBackwardInternal ::invoke(member, m, plen, piv.data(), ps, A.data(), as1, as0);
    }
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
