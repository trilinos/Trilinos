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

#ifndef KOKKOSSPARSE_SPMV_TEAM_SPEC_HPP_
#define KOKKOSSPARSE_SPMV_TEAM_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>
#include <KokkosSparse_spmv_team_impl.hpp>

namespace KokkosSparse {

template <typename MemberType>
struct TeamSpmv {
  template <typename ScalarType, typename ValuesViewType, typename IntView, typename xViewType, typename yViewType,
            int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha,
                                           const ValuesViewType& values, const IntView& row_ptr,
                                           const IntView& colIndices, const xViewType& x, const ScalarType beta,
                                           const yViewType& y) {
    return Impl::TeamSpmvInternal::invoke<MemberType, ScalarType, typename ValuesViewType::non_const_value_type,
                                          typename IntView::non_const_value_type, dobeta>(
        member, y.extent(0), alpha, values.data(), values.stride_0(), row_ptr.data(), row_ptr.stride_0(),
        colIndices.data(), colIndices.stride_0(), x.data(), x.stride_0(), beta, y.data(), y.stride_0());
  }
};

template <typename MemberType>
struct TeamVectorSpmv {
  template <typename ScalarType, typename ValuesViewType, typename IntView, typename xViewType, typename yViewType,
            int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha,
                                           const ValuesViewType& values, const IntView& row_ptr,
                                           const IntView& colIndices, const xViewType& x, const ScalarType beta,
                                           const yViewType& y) {
    return Impl::TeamVectorSpmvInternal::invoke<MemberType, ScalarType, typename ValuesViewType::non_const_value_type,
                                                typename IntView::non_const_value_type, dobeta>(
        member, y.extent(0), alpha, values.data(), values.stride_0(), row_ptr.data(), row_ptr.stride_0(),
        colIndices.data(), colIndices.stride_0(), x.data(), x.stride_0(), beta, y.data(), y.stride_0());
  }
};

}  // namespace KokkosSparse

#endif
