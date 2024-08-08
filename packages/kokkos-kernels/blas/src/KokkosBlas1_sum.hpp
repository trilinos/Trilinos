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

#ifndef KOKKOSBLAS1_SUM_HPP_
#define KOKKOSBLAS1_SUM_HPP_

#include <KokkosBlas1_sum_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief Return the sum of the vector x.
///
/// \tparam execution_space a Kokkos execution space where the kernel will run.
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param space [in] execution space instance where the kernel will run.
/// \param x [in] Input 1-D View.
///
/// \return The sum product result; a single value.
template <class execution_space, class XVector,
          typename std::enable_if<Kokkos::is_execution_space_v<execution_space>, int>::type = 0>
typename XVector::non_const_value_type sum(const execution_space& space, const XVector& x) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::sum: execution_space must be a valid Kokkos "
                "execution space");
  static_assert(Kokkos::is_view<XVector>::value, "KokkosBlas::sum: XVector must be a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible,
                "KokkosBlas::sum: XVector must be accessible from execution_space.");
  static_assert(XVector::rank == 1,
                "KokkosBlas::sum: "
                "Both Vector inputs must have rank 1.");

  using XVector_Internal = Kokkos::View<typename XVector::const_value_type*,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
                                        typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  using layout_t = typename XVector_Internal::array_layout;

  using RVector_Internal = Kokkos::View<typename XVector::non_const_value_type, layout_t, Kokkos::HostSpace,
                                        Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  typename XVector::non_const_value_type result;
  RVector_Internal R = RVector_Internal(&result, layout_t());
  XVector_Internal X = x;

  Impl::Sum<execution_space, RVector_Internal, XVector_Internal>::sum(space, R, X);
  space.fence();
  return result;
}

/// \brief Return the sum of the vector x.
///
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XVector.
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
///
/// \return The sum product result; a single value.
template <class XVector>
typename XVector::non_const_value_type sum(const XVector& x) {
  return sum(typename XVector::execution_space{}, x);
}

/// \brief R(j) = sum(X(i,j))
///
/// Replace each entry in R with the sumolute value (magnitude) of the
/// corresponding entry in X.
/// This function is non-blocking and thread-safe.
///
/// \tparam execution_space a Kokkos execution space where the kernel will run.
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
///
/// \param space [in] execution space instance where the kernel will run.
/// \param R [out] Output View (rank 0 or 1) containing the results.
/// \param X [in] Input View (rank 1 or 2).
template <class execution_space, class RV, class XMV>
void sum(const execution_space& space, const RV& R, const XMV& X,
         typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::sum: execution_space must be a valid Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::sum: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::sum: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::sum: XMV must be accessible from execution_space.");
  static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                "KokkosBlas::sum: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert(((RV::rank == 0) && (XMV::rank == 1)) || ((RV::rank == 1) && (XMV::rank == 2)),
                "KokkosBlas::sum: "
                "RV and XMV must either have rank 0 and 1 or rank 1 and 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(1) != R.extent(0)) {
    std::ostringstream os;
    os << "KokkosBlas::sum (MV): Dimensions of R and X do not match: "
       << "R: " << R.extent(0) << ", X: " << X.extent(0) << " x " << X.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using UnifiedXLayout  = typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using UnifiedRVLayout = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<RV, UnifiedXLayout>::array_layout;

  // Create unmanaged versions of the input Views.  RV and XMV may be
  // rank 1 or rank 2.
  typedef Kokkos::View<typename RV::non_const_data_type, UnifiedRVLayout, typename RV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      RV_Internal;
  typedef Kokkos::View<typename XMV::const_data_type, UnifiedXLayout, typename XMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XMV_Internal;

  RV_Internal R_internal  = R;
  XMV_Internal X_internal = X;

  Impl::Sum<execution_space, RV_Internal, XMV_Internal>::sum(space, R_internal, X_internal);
}

/// \brief R(j) = sum(X(i,j))
///
/// Replace each entry in R with the sumolute value (magnitude) of the
/// corresponding entry in X.
/// This function is non-blocking and thread-safe.
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XVM.
///
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
///
/// \param R [out] Output View (rank 0 or 1) containing the results.
/// \param X [in] Input View (rank 1 or 2).
template <class RV, class XMV>
void sum(const RV& R, const XMV& X, typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  sum(typename XMV::execution_space{}, R, X);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_SUM_HPP_
