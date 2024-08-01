/*
//@HEADER
// *****************************************************************************
//                        Adelus
//
// Copyright 2020 NTESS and the Adelus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
*/

#ifndef BLASWRAPPER_COPY_HPP_
#define BLASWRAPPER_COPY_HPP_

#include <stdexcept>

#include<BlasWrapper_copy_spec.hpp>
#include<KokkosKernels_helpers.hpp>

namespace BlasWrapper {

/// \brief Y(i,j) = X(i,j)
///
/// Copy each entry in X into the corresponding entry in Y.
///
/// \tparam YMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as YMV, and its entries must be assignable to
///   those of YMV.
template<class XMV, class YMV>
void
copy (const XMV& X, const YMV& Y)
{
  static_assert (Kokkos::is_view<XMV>::value, "BlasWrapper::copy: "
                 "X is not a Kokkos::View.");
  static_assert (Kokkos::is_view<YMV>::value, "BlasWrapper::copy: "
                 "Y is not a Kokkos::View.");
  static_assert (std::is_same<typename YMV::value_type,
                 typename YMV::non_const_value_type>::value,
                 "BlasWrapper::copy: Y is const.  "
                 "It must be nonconst, because it is an output argument "
                 "(we have to be able to write to its entries).");
  static_assert (int(XMV::rank) == int(YMV::rank), "BlasWrapper::copy: "
                 "X and Y must have the same rank.");
  static_assert (YMV::rank == 1 || YMV::rank == 2, "BlasWrapper::copy: "
                 "YMV and XMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != Y.extent(0) ||
      X.extent(1) != Y.extent(1)) {
    std::ostringstream os;
    os << "BlasWrapper::copy (MV): Dimensions of Y and X do not match: "
       << "Y: " << Y.extent(0) << " x " << Y.extent(1)
       << ", X: " << X.extent(0) << " x " << X.extent(1);
    throw std::runtime_error (os.str ());
  }

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.
  typedef Kokkos::View<
    typename std::conditional<
      YMV::rank == 1,
      typename YMV::non_const_value_type*,
      typename YMV::non_const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<YMV>::array_layout,
    typename YMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > YMV_Internal;
  typedef Kokkos::View<
    typename std::conditional<
      XMV::rank == 1,
      typename XMV::const_value_type*,
      typename XMV::const_value_type** >::type,
    typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout,
    typename XMV::device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged> > XMV_Internal;

  YMV_Internal Y_internal = Y;
  XMV_Internal X_internal = X;

  Impl::Copy<XMV_Internal, YMV_Internal>::copy (X_internal, Y_internal);
}
}

#endif // BLASWRAPPER_COPY_HPP_
