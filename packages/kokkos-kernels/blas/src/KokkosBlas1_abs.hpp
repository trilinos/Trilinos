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

#ifndef KOKKOSBLAS1_ABS_HPP_
#define KOKKOSBLAS1_ABS_HPP_

#include <KokkosBlas1_abs_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief R(i,j) = abs(X(i,j))
///
/// Non-blocking function to replace each entry in R with the absolute value
/// (magnitude) of the corresponding entry in X.
///
/// \tparam execution_space a Kokkos execution space to run the kernels on.
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
///
/// \param space [in] an execution_space instance where the kernel will run.
/// \param R [out] view of type RMV that contains the absolute value X on
/// output.
/// \param X [in] view of type XMV.
template <class execution_space, class RMV, class XMV>
void abs(const execution_space& space, const RMV& R, const XMV& X) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::abs: execution_space must be a valid Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view<RMV>::value,
                "KokkosBlas::abs: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename RMV::memory_space>::accessible,
                "KokkosBlas::abs: RMV must be accessible from execution space");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::abs: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::abs: XMV must be accessible from execution space");
  static_assert(std::is_same<typename RMV::value_type, typename RMV::non_const_value_type>::value,
                "KokkosBlas::abs: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert(int(RMV::rank) == int(XMV::rank),
                "KokkosBlas::abs: "
                "R and X must have the same rank.");
  static_assert(RMV::rank == 1 || RMV::rank == 2,
                "KokkosBlas::abs: "
                "RMV and XMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != R.extent(0) || X.extent(1) != R.extent(1)) {
    std::ostringstream os;
    os << "KokkosBlas::abs (MV): Dimensions of R and X do not match: "
       << "R: " << R.extent(0) << " x " << R.extent(1) << ", X: " << X.extent(0) << " x " << X.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.
  using RMV_Internal = Kokkos::View<typename std::conditional<RMV::rank == 1, typename RMV::non_const_value_type*,
                                                              typename RMV::non_const_value_type**>::type,
                                    typename KokkosKernels::Impl::GetUnifiedLayout<RMV>::array_layout,
                                    typename RMV::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using XMV_Internal = Kokkos::View<typename std::conditional<XMV::rank == 1, typename XMV::const_value_type*,
                                                              typename XMV::const_value_type**>::type,
                                    typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout,
                                    typename XMV::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  RMV_Internal R_internal = R;
  XMV_Internal X_internal = X;

  Impl::Abs<execution_space, RMV_Internal, XMV_Internal>::abs(space, R_internal, X_internal);
}

/// \brief R(i,j) = abs(X(i,j))
///
/// Non-blocking function to replace each entry in R with the absolute value
/// (magnitude) of the corresponding entry in X. The kernel is executed in the
/// default stream/queue associated with the execution space of RMV.
///
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
///
/// \param R [out] view of type RMV that contains the absolute value X on
/// output.
/// \param X [in] view of type XMV.
template <class RMV, class XMV>
void abs(const RMV& R, const XMV& X) {
  abs(typename RMV::execution_space{}, R, X);
}
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_ABS_HPP_
