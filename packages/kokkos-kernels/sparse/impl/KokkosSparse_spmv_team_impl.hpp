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

#ifndef KOKKOSSPARSE_SPMV_TEAM_IMPL_HPP_
#define KOKKOSSPARSE_SPMV_TEAM_IMPL_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosSparse {
namespace Impl {

struct TeamSpmvInternal {
  template <typename MemberType, typename ScalarType, typename ValueType, typename OrdinalType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OrdinalType numRows, const ScalarType alpha,
                                           const ValueType* KOKKOS_RESTRICT values, const OrdinalType valuess0,
                                           const OrdinalType* KOKKOS_RESTRICT row_ptr, const OrdinalType row_ptrs0,
                                           const OrdinalType* KOKKOS_RESTRICT colIndices,
                                           const OrdinalType colIndicess0, const ValueType* KOKKOS_RESTRICT x,
                                           const OrdinalType xs0, const ScalarType beta,
                                           /**/ ValueType* KOKKOS_RESTRICT y, const OrdinalType ys0);
};

struct TeamVectorSpmvInternal {
  template <typename MemberType, typename ScalarType, typename ValueType, typename OrdinalType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const OrdinalType numRows, const ScalarType alpha,
                                           const ValueType* KOKKOS_RESTRICT values, const OrdinalType valuess0,
                                           const OrdinalType* KOKKOS_RESTRICT row_ptr, const OrdinalType row_ptrs0,
                                           const OrdinalType* KOKKOS_RESTRICT colIndices,
                                           const OrdinalType colIndicess0, const ValueType* KOKKOS_RESTRICT x,
                                           const OrdinalType xs0, const ScalarType beta,
                                           /**/ ValueType* KOKKOS_RESTRICT y, const OrdinalType ys0);
};

template <typename MemberType, typename ScalarType, typename ValueType, typename OrdinalType, int dobeta>
KOKKOS_INLINE_FUNCTION int TeamSpmvInternal::invoke(
    const MemberType& member, const OrdinalType numRows, const ScalarType alpha,
    const ValueType* KOKKOS_RESTRICT values, const OrdinalType valuess0, const OrdinalType* KOKKOS_RESTRICT row_ptr,
    const OrdinalType row_ptrs0, const OrdinalType* KOKKOS_RESTRICT colIndices, const OrdinalType colIndicess0,
    const ValueType* KOKKOS_RESTRICT x, const OrdinalType xs0, const ScalarType beta,
    /**/ ValueType* KOKKOS_RESTRICT y, const OrdinalType ys0) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numRows), [&](const OrdinalType& iRow) {
    const OrdinalType rowLength = row_ptr[(iRow + 1) * row_ptrs0] - row_ptr[iRow * row_ptrs0];
    ValueType sum               = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (OrdinalType iEntry = 0; iEntry < rowLength; ++iEntry) {
      sum += values[(row_ptr[iRow * row_ptrs0] + iEntry) * valuess0] *
             x[colIndices[(row_ptr[iRow * row_ptrs0] + iEntry) * colIndicess0] * xs0];
    }

    sum *= alpha;

    if (dobeta == 0) {
      y[iRow * ys0] = sum;
    } else {
      y[iRow * ys0] = beta * y[iRow * ys0] + sum;
    }
  });
  return 0;
}

template <typename MemberType, typename ScalarType, typename ValueType, typename OrdinalType, int dobeta>
KOKKOS_INLINE_FUNCTION int TeamVectorSpmvInternal::invoke(
    const MemberType& member, const OrdinalType numRows, const ScalarType alpha,
    const ValueType* KOKKOS_RESTRICT values, const OrdinalType valuess0, const OrdinalType* KOKKOS_RESTRICT row_ptr,
    const OrdinalType row_ptrs0, const OrdinalType* KOKKOS_RESTRICT colIndices, const OrdinalType colIndicess0,
    const ValueType* KOKKOS_RESTRICT x, const OrdinalType xs0, const ScalarType beta,
    /**/ ValueType* KOKKOS_RESTRICT y, const OrdinalType ys0) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, numRows), [&](const OrdinalType& iRow) {
    const OrdinalType rowLength = row_ptr[(iRow + 1) * row_ptrs0] - row_ptr[iRow * row_ptrs0];

    ValueType sum = 0;
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(member, rowLength),
        [&](const OrdinalType& iEntry, ValueType& val) {
          val += values[(row_ptr[iRow * row_ptrs0] + iEntry) * valuess0] *
                 x[colIndices[(row_ptr[iRow * row_ptrs0] + iEntry) * colIndicess0] * xs0];
        },
        sum);

    sum *= alpha;

    if (dobeta == 0) {
      y[iRow * ys0] = sum;
    } else {
      y[iRow * ys0] = beta * y[iRow * ys0] + sum;
    }
  });
  return 0;
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif
