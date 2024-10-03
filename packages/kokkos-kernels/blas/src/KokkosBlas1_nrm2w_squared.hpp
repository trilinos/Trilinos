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

#ifndef KOKKOSBLAS1_NRM2W_SQUARED_HPP_
#define KOKKOSBLAS1_NRM2W_SQUARED_HPP_

#include <KokkosBlas1_nrm2w_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief Return the nrm2w of the vector x.
///
/// \tparam execution_space a Kokkos execution space where the computation
///         will be launched
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param space [in] an execution space instance that may specify
/// \param x [in] Input 1-D View.
/// \param w [in] Input weights (1-D View).
///
/// \return The nrm2w product result; a single value.
template <class execution_space, class XVector,
          typename std::enable_if<Kokkos::is_execution_space<execution_space>::value, int>::type = 0>
typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type nrm2w_squared(
    const execution_space& space, const XVector& x, const XVector& w) {
  static_assert(Kokkos::is_execution_space<execution_space>::value,
                "KokkosBlas::nrm2w_squared: execution_space must be a valid "
                "Kokkos execution space.");
  static_assert(Kokkos::is_view<XVector>::value, "KokkosBlas::nrm2w_squared: XVector must be a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XVector::memory_space>::accessible,
                "KokkosBlas::nrm2w_squared: XVector must be accessible from "
                "execution_space.");
  static_assert(XVector::rank == 1,
                "KokkosBlas::nrm2w_squared: "
                "Both Vector inputs must have rank 1.");

  using mag_type = typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type;

  using XVector_Internal = Kokkos::View<typename XVector::const_value_type*,
                                        typename KokkosKernels::Impl::GetUnifiedLayout<XVector>::array_layout,
                                        typename XVector::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  using layout_t = typename XVector_Internal::array_layout;

  using RVector_Internal =
      Kokkos::View<mag_type, layout_t, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  mag_type result;
  RVector_Internal R = RVector_Internal(&result, layout_t());
  XVector_Internal X = x;
  XVector_Internal W = w;

  Impl::Nrm2w<execution_space, RVector_Internal, XVector_Internal>::nrm2w(space, R, X, W, false);
  space.fence();
  return result;
}

/// \brief Return the nrm2w of the vector x.
///
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XVector.
///
/// \tparam XVector Type of the first vector x; a 1-D Kokkos::View.
///
/// \param x [in] Input 1-D View.
/// \param w [in] Input weights (1-D View).
///
/// \return The nrm2w product result; a single value.
template <class XVector>
typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::mag_type nrm2w_squared(
    const XVector& x, const XVector& w) {
  return nrm2w_squared(typename XVector::execution_space(), x, w);
}

/// \brief R(i,j) = nrm2w(X(i,j))
///
/// Replace each entry in R with the nrm2wolute value (magnitude) of the
/// corresponding entry in X.
/// This function is non-blocking and thread-safe.
///
/// \tparam execution_space a Kokkos execution space where the kernel will run.
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
///
/// \param space [in] the execution space instance, possibly containing a
/// stream/queue where the kernel will be executed.
/// \param R [out] Output View containing results (rank 0 or 1).
/// \param X [in] Input View (rank 1 or 2).
/// \param W [in] Input View (rank 1 or 2).
template <class execution_space, class RV, class XMV>
void nrm2w_squared(const execution_space& space, const RV& R, const XMV& X, const XMV& W,
                   typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  static_assert(Kokkos::is_execution_space<execution_space>::value,
                "KokkosBlas::nrm2w_squared: execution_space must be a valid "
                "Kokkos execution space.");
  static_assert(Kokkos::is_view<RV>::value,
                "KokkosBlas::nrm2w_squared: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::nrm2w_squared: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::nrm2w_squared: XVector must be accessible from "
                "execution_space.");
  static_assert(std::is_same<typename RV::value_type, typename RV::non_const_value_type>::value,
                "KokkosBlas::nrm2w_squared: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert(((RV::rank == 0) && (XMV::rank == 1)) || ((RV::rank == 1) && (XMV::rank == 2)),
                "KokkosBlas::nrm2w_squared: "
                "RV and XMV must either have rank 0 and 1 or rank 1 and 2.");
  using mag_type = typename Kokkos::Details::InnerProductSpaceTraits<typename XMV::non_const_value_type>::mag_type;
  static_assert(std::is_same<typename RV::value_type, mag_type>::value,
                "KokkosBlas::nrm2w: R must have the magnitude type of"
                "the xvectors value_type it is an output argument "
                "(we have to be able to write to its entries).");

  // Check compatibility of dimensions at run time.
  if (X.extent(1) != R.extent(0)) {
    std::ostringstream os;
    os << "KokkosBlas::nrm2w (MV): Dimensions of R and X do not match: "
       << "R: " << R.extent(0) << ", X: " << X.extent(0) << " x " << X.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using UnifiedXLayout  = typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using UnifiedRVLayout = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<RV, UnifiedXLayout>::array_layout;

  // Create unmanaged versions of the input Views.  RV and XMV may be
  // rank 1 or rank 2.
  using RV_Internal  = Kokkos::View<typename RV::non_const_data_type, UnifiedRVLayout, typename RV::device_type,
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using XMV_Internal = Kokkos::View<typename XMV::const_data_type, UnifiedXLayout, typename XMV::device_type,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  RV_Internal R_internal  = R;
  XMV_Internal X_internal = X;
  XMV_Internal W_internal = W;

  Impl::Nrm2w<execution_space, RV_Internal, XMV_Internal>::nrm2w(space, R_internal, X_internal, W_internal, false);
}

/// \brief R(i,j) = nrm2w(X(i,j))
///
/// Replace each entry in R with the nrm2wolute value (magnitude) of the
/// corresponding entry in X.
/// This function is non-blocking and thread-safe.
/// The kernel is executed in the default stream/queue associated
/// with the execution space of XMV.
///
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
///
/// \param R [out] Output View containing results (rank 0 or 1).
/// \param X [in] Input View (rank 1 or 2).
/// \param W [in] Input View (rank 1 or 2).
template <class RV, class XMV>
void nrm2w_squared(const RV& R, const XMV& X, const XMV& W,
                   typename std::enable_if<Kokkos::is_view<RV>::value, int>::type = 0) {
  nrm2w_squared(typename XMV::execution_space{}, R, X, W);
}
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_NRM2W_HPP_
