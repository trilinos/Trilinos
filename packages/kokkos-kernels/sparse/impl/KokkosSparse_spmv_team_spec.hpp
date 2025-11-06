// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSSPARSE_SPMV_TEAM_SPEC_HPP_
#define KOKKOSSPARSE_SPMV_TEAM_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_ArithTraits.hpp>
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
        member, y.extent(0), alpha, values.data(), values.stride(0), row_ptr.data(), row_ptr.stride(0),
        colIndices.data(), colIndices.stride(0), x.data(), x.stride(0), beta, y.data(), y.stride(0));
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
        member, y.extent(0), alpha, values.data(), values.stride(0), row_ptr.data(), row_ptr.stride(0),
        colIndices.data(), colIndices.stride(0), x.data(), x.stride(0), beta, y.data(), y.stride(0));
  }
};

}  // namespace KokkosSparse

#endif
