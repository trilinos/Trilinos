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

#ifndef KOKKOSBLAS1_AXPBY_HPP_
#define KOKKOSBLAS1_AXPBY_HPP_

#include <KokkosBlas1_axpby_spec.hpp>
#include <KokkosBlas_serial_axpy.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

// axpby() accepts both scalar coefficients a and b, and vector
// coefficients (apply one for each column of the input multivectors).
// This traits class helps axpby() select the correct specialization
// of AV and BV (the type of a resp. b) for invoking the
// implementation.

namespace KokkosBlas {

template <class AV, class XMV, class BV, class YMV>
void axpby(const AV& a, const XMV& X, const BV& b, const YMV& Y) {
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::axpby: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::is_view<YMV>::value,
                "KokkosBlas::axpby: "
                "Y is not a Kokkos::View.");
  static_assert(std::is_same<typename YMV::value_type,
                             typename YMV::non_const_value_type>::value,
                "KokkosBlas::axpby: Y is const.  It must be nonconst, "
                "because it is an output argument "
                "(we must be able to write to its entries).");
  static_assert(int(YMV::Rank) == int(XMV::Rank),
                "KokkosBlas::axpby: "
                "X and Y must have the same rank.");
  static_assert(YMV::Rank == 1 || YMV::Rank == 2,
                "KokkosBlas::axpby: "
                "XMV and YMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    std::ostringstream os;
    os << "KokkosBlas::axpby: Dimensions of X and Y do not match: "
       << "X: " << X.extent(0) << " x " << X.extent(1) << ", Y: " << Y.extent(0)
       << " x " << Y.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using UnifiedXLayout =
      typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using UnifiedYLayout =
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<
          YMV, UnifiedXLayout>::array_layout;

  // Create unmanaged versions of the input Views.  XMV and YMV may be
  // rank 1 or rank 2.  AV and BV may be either rank-1 Views, or
  // scalar values.
  typedef Kokkos::View<typename XMV::const_data_type, UnifiedXLayout,
                       typename XMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XMV_Internal;
  typedef Kokkos::View<typename YMV::non_const_data_type, UnifiedYLayout,
                       typename YMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      YMV_Internal;
  typedef typename KokkosKernels::Impl::GetUnifiedScalarViewType<
      AV, XMV_Internal, true>::type AV_Internal;
  typedef typename KokkosKernels::Impl::GetUnifiedScalarViewType<
      BV, YMV_Internal, true>::type BV_Internal;

  AV_Internal a_internal  = a;
  XMV_Internal X_internal = X;
  BV_Internal b_internal  = b;
  YMV_Internal Y_internal = Y;

  Impl::Axpby<AV_Internal, XMV_Internal, BV_Internal, YMV_Internal>::axpby(
      a_internal, X_internal, b_internal, Y_internal);
}

template <class AV, class XMV, class YMV>
void axpy(const AV& a, const XMV& X, const YMV& Y) {
  axpby(a, X,
        Kokkos::Details::ArithTraits<typename YMV::non_const_value_type>::one(),
        Y);
}

///
/// Serial axpy on device
///
template <class scalar_type, class XMV, class YMV>
KOKKOS_FUNCTION void serial_axpy(const scalar_type alpha, const XMV X, YMV Y) {
#if (KOKKOSKERNELS_DEBUG_LEVEL > 0)
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::serial_axpy: XMV is not a Kokkos::View");
  static_assert(Kokkos::is_view<YMV>::value,
                "KokkosBlas::serial_axpy: YMV is not a Kokkos::View");
  static_assert(XMV::Rank == 1 || XMV::Rank == 2,
                "KokkosBlas::serial_axpy: XMV must have rank 1 or 2.");
  static_assert(
      XMV::Rank == YMV::Rank,
      "KokkosBlas::serial_axpy: XMV and YMV must have the same rank.");

  if (X.extent(0) != Y.extent(0) || X.extent(1) != Y.extent(1)) {
    Kokkos::abort("KokkosBlas::serial_axpy: X and Y dimensions do not match");
  }
#endif  // KOKKOSKERNELS_DEBUG_LEVEL

  return Impl::serial_axpy_mv(X.extent(0), X.extent(1), alpha, X.data(),
                              Y.data(), X.stride_0(), X.stride_1(),
                              Y.stride_0(), Y.stride_1());
}

}  // namespace KokkosBlas

#endif
