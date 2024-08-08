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
#ifndef __KOKKOSBATCHED_FIND_AMAX_INTERNAL_HPP__
#define __KOKKOSBATCHED_FIND_AMAX_INTERNAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// =====================
struct SerialFindAmaxInternal {
  template <typename ValueType, typename IntType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ValueType *KOKKOS_RESTRICT A, const int as0,
                                           /**/ IntType *KOKKOS_RESTRICT idx) {
    ValueType max_val(A[0]);
    IntType val_loc(0);
    for (int i = 1; i < m; ++i) {
      const int idx_a = i * as0;
      if (A[idx_a] > max_val) {
        max_val = A[idx_a];
        val_loc = i;
      }
    }
    *idx = val_loc;
    return 0;
  }
};

///
/// TeamVector Internal Impl
/// ========================
struct TeamVectorFindAmaxInternal {
  template <typename MemberType, typename ValueType, typename IntType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType &member, const int m, const ValueType *KOKKOS_RESTRICT A,
                                           const int as0,
                                           /**/ IntType *KOKKOS_RESTRICT idx) {
    if (m > 0) {
      using reducer_value_type = typename Kokkos::MaxLoc<ValueType, IntType>::value_type;
      reducer_value_type value{};
      Kokkos::MaxLoc<ValueType, IntType> reducer_value(value);
      Kokkos::parallel_reduce(
          Kokkos::TeamVectorRange(member, m),
          [&](const int &i, reducer_value_type &update) {
            const int idx_a = i * as0;
            if (A[idx_a] > update.val) {
              update.val = A[idx_a];
              update.loc = i;
            }
          },
          reducer_value);
      Kokkos::single(Kokkos::PerTeam(member), [&]() { *idx = value.loc; });
    } else {
      Kokkos::single(Kokkos::PerTeam(member), [&]() { *idx = 0; });
    }
    return 0;
  }
};

}  // namespace KokkosBatched

#endif
