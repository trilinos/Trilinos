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

#ifndef KOKKOSBLAS1_IAMAX_HPP_
#define KOKKOSBLAS1_IAMAX_HPP_

#include <KokkosBlas1_iamax_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief Return the (smallest) index of the element of the maximum magnitude
/// of the vector x.
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
///
/// \return The (smallest) index of the element of the maximum magnitude; a
/// single value.
///         Note: Returned index is 1-based for compatibility with Fortran.
template <class XVector>
typename XVector::size_type iamax(const XVector& x) {
  static_assert(Kokkos::is_view<XVector>::value,
                "KokkosBlas::iamax: XVector must be a Kokkos::View.");
  static_assert(XVector::rank == 1,
                "KokkosBlas::iamax: "
                "Both Vector inputs must have rank 1.");

  typedef typename XVector::size_type index_type;

  typedef Kokkos::View<
      typename XVector::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
      typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XVector_Internal;

  typedef Kokkos::View<index_type, typename XVector_Internal::array_layout,
                       Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      RVector_Internal;

  index_type result;
  RVector_Internal R = RVector_Internal(&result);
  XVector_Internal X = x;

  Impl::Iamax<RVector_Internal, XVector_Internal>::iamax(R, X);
  Kokkos::fence();
  return result;
}

/// \brief R(j) = iamax(X(i,j))
///
/// Replace each entry in R with the (smallest) index of the element of the
/// maximum magnitude of the corresponding entry in X.
///
/// \tparam RMV 0-D or 1-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.
///
/// Note for TPL cuBLAS: When TPL cuBLAS iamax is used and returns result to a
/// view, RMV must be 0-D view and XMV must be 1-D view.
template <class RV, class XMV>
void iamax(const RV& R, const XMV& X,
           typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::iamax: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::iamax: "
                "X is not a Kokkos::View.");
  static_assert(std::is_same<typename RV::value_type,
                             typename RV::non_const_value_type>::value,
                "KokkosBlas::iamax: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert(((RV::rank == 0) && (XMV::rank == 1)) ||
                    ((RV::rank == 1) && (XMV::rank == 2)),
                "KokkosBlas::iamax: "
                "RV and XMV must either have rank 0 and 1 or rank 1 and 2.");

  typedef typename XMV::size_type index_type;
  static_assert(std::is_same<typename RV::value_type, index_type>::value,
                "KokkosBlas::iamax: R must have the type of"
                "the Xvectors size_type it is an output argument "
                "(we have to be able to write to its entries).");

  // Check compatibility of dimensions at run time.
  if (X.extent(1) != R.extent(0)) {
    std::ostringstream os;
    os << "KokkosBlas::iamax (MV): Dimensions of R and X do not match: "
       << "R: " << R.extent(0) << ", X: " << X.extent(0) << " x "
       << X.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using UnifiedXLayout =
      typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using UnifiedRVLayout =
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<
          RV, UnifiedXLayout>::array_layout;

  // Create unmanaged versions of the input Views.  RV may be rank 0 or rank 1.
  // XMV may be rank 1 or rank 2.
  typedef Kokkos::View<typename std::conditional<
                           RV::rank == 0, typename RV::non_const_value_type,
                           typename RV::non_const_value_type*>::type,
                       UnifiedRVLayout,
                       typename std::conditional<
                           std::is_same<typename RV::device_type::memory_space,
                                        Kokkos::HostSpace>::value,
                           Kokkos::HostSpace, typename RV::device_type>::type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      RV_Internal;
  typedef Kokkos::View<
      typename std::conditional<XMV::rank == 1, typename XMV::const_value_type*,
                                typename XMV::const_value_type**>::type,
      UnifiedXLayout, typename XMV::device_type,
      Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XMV_Internal;

  RV_Internal R_internal  = R;
  XMV_Internal X_internal = X;

  Impl::Iamax<RV_Internal, XMV_Internal>::iamax(R_internal, X_internal);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_IAMAX_HPP_
