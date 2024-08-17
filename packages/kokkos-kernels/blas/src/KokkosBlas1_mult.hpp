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

#ifndef KOKKOSBLAS1_MULT_HPP_
#define KOKKOSBLAS1_MULT_HPP_

#include <KokkosBlas1_mult_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief Element wise multiplication of two vectors:
///        Y[i] = gamma * Y[i] + alpha * A[i] * X[i]
///
/// This function is non-blocking and thread-safe
///
/// \tparam execution_type a Kokkos execution space type.
/// \tparam YMV Type of the first vector Y; a 1-D or 2-D Kokkos::View.
/// \tparam AV  Type of the second vector A; a 1-D Kokkos::View.
/// \tparam XMV Type of the third vector X; a 1-D or 2-D Kokkos::View.
///
/// \param space [in] An instance of execution_space on which the kernel
///                   will run (it may specify an execution stream/queue).
/// \param gamma [in] The scalar to apply to Y.
/// \param Y [in/out] The Y vector.
/// \param alpha [in] The scalar to apply to A.
/// \param A [in]     The vector to apply to X.
/// \param X [in]     The X vector.
template <class execution_space, class YMV, class AV, class XMV>
void mult(const execution_space& space, typename YMV::const_value_type& gamma, const YMV& Y,
          typename AV::const_value_type& alpha, const AV& A, const XMV& X) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::mult: execution_space must be a valid Kokkos "
                "execution space.");
  static_assert(Kokkos::is_view<YMV>::value,
                "KokkosBlas::mult: "
                "Y is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename YMV::memory_space>::accessible,
                "KokkosBlas::mult: YMV must be accessible from execution_space.");
  static_assert(Kokkos::is_view<AV>::value,
                "KokkosBlas::mult: "
                "A is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename AV::memory_space>::accessible,
                "KokkosBlas::mult: AV must be accessible from execution_space.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::mult: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::mult: AV must be accessible from execution_space.");
  static_assert(std::is_same<typename YMV::value_type, typename YMV::non_const_value_type>::value,
                "KokkosBlas::mult: Y is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert((XMV::rank == 1 && YMV::rank == 1) || (XMV::rank == 2 && YMV::rank == 2),
                "KokkosBlas::mult: Y and X must be either both rank 1, "
                "or both rank 2.");
  static_assert(AV::rank == 1, "KokkosBlas::mult: A must have rank 1.");

  // Check compatibility of dimensions at run time.
  if (Y.extent(0) != A.extent(0) || Y.extent(0) != X.extent(0) || Y.extent(1) != X.extent(1)) {
    std::ostringstream os;
    os << "KokkosBlas::mult: Dimensions do not match: "
       << "Y: " << Y.extent(0) << " x " << Y.extent(1) << ", A: " << A.extent(0) << " x " << A.extent(0)
       << ", X: " << X.extent(0) << " x " << X.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using YUnifiedLayout = typename KokkosKernels::Impl::GetUnifiedLayout<YMV>::array_layout;
  using AUnifiedLayout = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<AV, YUnifiedLayout>::array_layout;
  using XUnifiedLayout = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<XMV, YUnifiedLayout>::array_layout;

  // Create unmanaged versions of the input Views.
  typedef Kokkos::View<typename YMV::non_const_data_type, YUnifiedLayout, typename YMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      YMV_Internal;
  typedef Kokkos::View<typename AV::const_value_type*, AUnifiedLayout, typename AV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      AV_Internal;
  typedef Kokkos::View<typename XMV::const_data_type, XUnifiedLayout, typename XMV::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      XMV_Internal;

  YMV_Internal Y_internal = Y;
  AV_Internal A_internal  = A;
  XMV_Internal X_internal = X;

  Impl::Mult<execution_space, YMV_Internal, AV_Internal, XMV_Internal>::mult(space, gamma, Y_internal, alpha,
                                                                             A_internal, X_internal);
}

/// \brief Element wise multiplication of two vectors:
///        Y[i] = gamma * Y[i] + alpha * A[i] * X[i]
///
/// This function is non-blocking and thread-safe
/// The kernel is executed in the default stream/queue
/// associated with the execution space of YMV.
///
/// \tparam YMV Type of the first vector Y; a 1-D or 2-D Kokkos::View.
/// \tparam AV  Type of the second vector A; a 1-D Kokkos::View.
/// \tparam XMV Type of the third vector X; a 1-D or 2-D Kokkos::View.
///
/// \param gamma [in] The scalar to apply to Y.
/// \param Y [in/out] The Y vector.
/// \param alpha [in] The scalar to apply to A.
/// \param A [in]     The vector to apply to X.
/// \param X [in]     The X vector.
template <class YMV, class AV, class XMV>
void mult(typename YMV::const_value_type& gamma, const YMV& Y, typename AV::const_value_type& alpha, const AV& A,
          const XMV& X) {
  mult(typename YMV::execution_space{}, gamma, Y, alpha, A, X);
}

}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_MULT_HPP_
