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
#ifndef __KOKKOSBATCHED_SPMV_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_SPMV_SERIAL_IMPL_HPP__

/// \author Kim Liegeois (knliege@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================
struct SerialSpmvInternal {
  template <typename ScalarType, typename ValueType, typename OrdinalType, typename layout, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const OrdinalType numMatrices, const OrdinalType numRows, const ScalarType* KOKKOS_RESTRICT alpha,
      const OrdinalType alphas0, const ValueType* KOKKOS_RESTRICT values, const OrdinalType valuess0,
      const OrdinalType valuess1, const OrdinalType* KOKKOS_RESTRICT row_ptr, const OrdinalType row_ptrs0,
      const OrdinalType* KOKKOS_RESTRICT colIndices, const OrdinalType colIndicess0, const ValueType* KOKKOS_RESTRICT X,
      const OrdinalType xs0, const OrdinalType xs1, const ScalarType* KOKKOS_RESTRICT beta, const OrdinalType betas0,
      /**/ ValueType* KOKKOS_RESTRICT Y, const OrdinalType ys0, const OrdinalType ys1) {
    for (OrdinalType iMatrix = 0; iMatrix < numMatrices; ++iMatrix) {
      for (OrdinalType iRow = 0; iRow < numRows; ++iRow) {
        const OrdinalType rowLength = row_ptr[(iRow + 1) * row_ptrs0] - row_ptr[iRow * row_ptrs0];
        ValueType sum               = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
        for (OrdinalType iEntry = 0; iEntry < rowLength; ++iEntry) {
          sum += values[iMatrix * valuess0 + (row_ptr[iRow * row_ptrs0] + iEntry) * valuess1] *
                 X[iMatrix * xs0 + colIndices[(row_ptr[iRow * row_ptrs0] + iEntry) * colIndicess0] * xs1];
        }

        sum *= alpha[iMatrix * alphas0];

        if (dobeta == 0) {
          Y[iMatrix * ys0 + iRow * ys1] = sum;
        } else {
          Y[iMatrix * ys0 + iRow * ys1] = beta[iMatrix * betas0] * Y[iMatrix * ys0 + iRow * ys1] + sum;
        }
      }
    }

    return 0;
  }

  template <typename ScalarType, typename ValueType, typename OrdinalType, typename layout, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const OrdinalType numMatrices, const OrdinalType numRows,
                                           const ScalarType alpha, const ValueType* KOKKOS_RESTRICT values,
                                           const OrdinalType valuess0, const OrdinalType valuess1,
                                           const OrdinalType* KOKKOS_RESTRICT row_ptr, const OrdinalType row_ptrs0,
                                           const OrdinalType* KOKKOS_RESTRICT colIndices,
                                           const OrdinalType colIndicess0, const ValueType* KOKKOS_RESTRICT X,
                                           const OrdinalType xs0, const OrdinalType xs1, const ScalarType beta,
                                           /**/ ValueType* KOKKOS_RESTRICT Y, const OrdinalType ys0,
                                           const OrdinalType ys1) {
    for (OrdinalType iMatrix = 0; iMatrix < numMatrices; ++iMatrix) {
      for (OrdinalType iRow = 0; iRow < numRows; ++iRow) {
        const OrdinalType rowLength = row_ptr[(iRow + 1) * row_ptrs0] - row_ptr[iRow * row_ptrs0];
        ValueType sum               = 0;
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
        for (OrdinalType iEntry = 0; iEntry < rowLength; ++iEntry) {
          sum += values[iMatrix * valuess0 + (row_ptr[iRow * row_ptrs0] + iEntry) * valuess1] *
                 X[iMatrix * xs0 + colIndices[(row_ptr[iRow * row_ptrs0] + iEntry) * colIndicess0] * xs1];
        }

        sum *= alpha;

        if (dobeta == 0) {
          Y[iMatrix * ys0 + iRow * ys1] = sum;
        } else {
          Y[iMatrix * ys0 + iRow * ys1] = beta * Y[iMatrix * ys0 + iRow * ys1] + sum;
        }
      }
    }

    return 0;
  }
};

template <>
struct SerialSpmv<Trans::NoTranspose> {
  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, typename alphaViewType,
            typename betaViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(const alphaViewType& alpha, const ValuesViewType& values,
                                           const IntView& row_ptr, const IntView& colIndices, const xViewType& X,
                                           const betaViewType& beta, const yViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<ValuesViewType>::value, "KokkosBatched::spmv: ValuesViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<IntView>::value, "KokkosBatched::spmv: IntView is not a Kokkos::View.");
    static_assert(Kokkos::is_view<xViewType>::value, "KokkosBatched::spmv: xViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<yViewType>::value, "KokkosBatched::spmv: yViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<alphaViewType>::value, "KokkosBatched::spmv: alphaViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<betaViewType>::value, "KokkosBatched::spmv: betaViewType is not a Kokkos::View.");

    static_assert(ValuesViewType::rank == 2, "KokkosBatched::spmv: ValuesViewType must have rank 2.");
    static_assert(IntView::rank == 1, "KokkosBatched::spmv: IntView must have rank 2.");
    static_assert(xViewType::rank == 2, "KokkosBatched::spmv: xViewType must have rank 2.");
    static_assert(yViewType::rank == 2, "KokkosBatched::spmv: yViewType must have rank 2.");
    static_assert(alphaViewType::rank == 1, "KokkosBatched::spmv: alphaViewType must have rank 1.");
    static_assert(betaViewType::rank == 1, "KokkosBatched::spmv: betaViewType must have rank 1.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::spmv: Dimensions of X and Y do not match: X: %d x "
          "%d, Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(0) != alpha.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::spmv: First dimension of X and alpha do not match: "
          "X: %d x %d, alpha: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)alpha.extent(0));
      return 1;
    }
    if (X.extent(0) != beta.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::spmv: First dimension of X and beta do not match: X: "
          "%d x %d, beta: %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)beta.extent(0));
      return 1;
    }
    if (X.extent(0) != values.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::spmv: First dimension of X and the first dimension "
          "of values do not match: X: %d x %d, values: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)values.extent(0), (int)values.extent(1));
      return 1;
    }
    if (colIndices.extent(0) != values.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::spmv: Dimension of colIndices and the second "
          "dimension of values do not match: colIndices: %d , values: %d x "
          "%d\n",
          (int)colIndices.extent(0), (int)values.extent(0), (int)values.extent(1));
      return 1;
    }
    if (row_ptr.extent(0) - 1 != X.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::spmv: Dimension of row_ptr and the second dimension "
          "of X do not match: colIndices (-1): %d , values: %d x %d\n",
          (int)row_ptr.extent(0) - 1, (int)X.extent(0), (int)X.extent(1));
      return 1;
    }
#endif

    return SerialSpmvInternal::template invoke<
        typename alphaViewType::non_const_value_type, typename ValuesViewType::non_const_value_type,
        typename IntView::non_const_value_type, typename ValuesViewType::array_layout, dobeta>(
        X.extent(0), X.extent(1), alpha.data(), alpha.stride_0(), values.data(), values.stride_0(), values.stride_1(),
        row_ptr.data(), row_ptr.stride_0(), colIndices.data(), colIndices.stride_0(), X.data(), X.stride_0(),
        X.stride_1(), beta.data(), beta.stride_0(), Y.data(), Y.stride_0(), Y.stride_1());
  }

  template <typename ValuesViewType, typename IntView, typename xViewType, typename yViewType, int dobeta>
  KOKKOS_INLINE_FUNCTION static int invoke(
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type& alpha,
      const ValuesViewType& values, const IntView& row_ptr, const IntView& colIndices, const xViewType& X,
      const typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type& beta,
      const yViewType& Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
    static_assert(Kokkos::is_view<ValuesViewType>::value, "KokkosBatched::spmv: ValuesViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<IntView>::value, "KokkosBatched::spmv: IntView is not a Kokkos::View.");
    static_assert(Kokkos::is_view<xViewType>::value, "KokkosBatched::spmv: xViewType is not a Kokkos::View.");
    static_assert(Kokkos::is_view<yViewType>::value, "KokkosBatched::spmv: yViewType is not a Kokkos::View.");

    static_assert(ValuesViewType::rank == 2, "KokkosBatched::spmv: ValuesViewType must have rank 2.");
    static_assert(IntView::rank == 1, "KokkosBatched::spmv: IntView must have rank 2.");
    static_assert(xViewType::rank == 2, "KokkosBatched::spmv: xViewType must have rank 2.");
    static_assert(yViewType::rank == 2, "KokkosBatched::spmv: yViewType must have rank 2.");

    // Check compatibility of dimensions at run time.
    if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::spmv: Dimensions of X and Y do not match: X: %d x "
          "%d, Y: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)Y.extent(0), (int)Y.extent(1));
      return 1;
    }
    if (X.extent(0) != values.extent(0)) {
      Kokkos::printf(
          "KokkosBatched::spmv: First dimension of X and the first dimension "
          "of values do not match: X: %d x %d, values: %d x %d\n",
          (int)X.extent(0), (int)X.extent(1), (int)values.extent(0), (int)values.extent(1));
      return 1;
    }
    if (colIndices.extent(0) != values.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::spmv: Dimension of colIndices and the second "
          "dimension of values do not match: colIndices: %d , values: %d x "
          "%d\n",
          (int)colIndices.extent(0), (int)values.extent(0), (int)values.extent(1));
      return 1;
    }
    if (row_ptr.extent(0) - 1 != X.extent(1)) {
      Kokkos::printf(
          "KokkosBatched::spmv: Dimension of row_ptr and the second dimension "
          "of X do not match: colIndices (-1): %d , values: %d x %d\n",
          (int)row_ptr.extent(0) - 1, (int)X.extent(0), (int)X.extent(1));
      return 1;
    }
#endif

    return SerialSpmvInternal::template invoke<
        typename Kokkos::ArithTraits<typename ValuesViewType::non_const_value_type>::mag_type,
        typename ValuesViewType::non_const_value_type, typename IntView::non_const_value_type,
        typename ValuesViewType::array_layout, dobeta>(
        X.extent(0), X.extent(1), alpha, values.data(), values.stride_0(), values.stride_1(), row_ptr.data(),
        row_ptr.stride_0(), colIndices.data(), colIndices.stride_0(), X.data(), X.stride_0(), X.stride_1(), beta,
        Y.data(), Y.stride_0(), Y.stride_1());
  }
};

}  // namespace KokkosBatched

#endif
