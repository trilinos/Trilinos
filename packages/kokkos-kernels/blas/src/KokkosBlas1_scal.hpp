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

#ifndef KOKKOSBLAS1_SCAL_HPP_
#define KOKKOSBLAS1_SCAL_HPP_

#include <KokkosBlas1_scal_spec.hpp>
#include <KokkosBlas1_serial_scal_impl.hpp>
#include <KokkosBlas1_team_scal_impl.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

///
/// General/Host Scale
///

namespace KokkosBlas {

/// \brief Computes R := alpha*X
///
/// This function is non-blocking and thread-safe
///
/// \tparam execution_space a Kokkos execution space where the kernel will run.
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization. It must have
///   the same rank as RMV.
/// \tparam AV 1-D or 2-D Kokkos::View specialization.
///
/// \param space [in] the execution space instance on which the kernel will run.
/// \param R [in/out] view of type RMV in which the results will be stored.
/// \param a [in] view of type AV, scaling parameter for X.
/// \param X [in] input view of type XMV.
template <class execution_space, class RMV, class AV, class XMV>
void scal(const execution_space& space, const RMV& R, const AV& a, const XMV& X) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::scal: execution_space must be a valid Kokkos "
                "execution space");
  static_assert(Kokkos::is_view<RMV>::value,
                "KokkosBlas::scal: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename RMV::memory_space>::accessible,
                "KokkosBlas::scal: RMV must be accessible from execution_space.");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::scal: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::scal: XMV must be accessible from execution_space");
  static_assert(Kokkos::SpaceAccessibility<typename RMV::memory_space, typename XMV::memory_space>::assignable,
                "KokkosBlas::scal: XMV must be assignable to RMV");
  static_assert(std::is_same<typename RMV::value_type, typename RMV::non_const_value_type>::value,
                "KokkosBlas::scal: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert((int)RMV::rank == (int)XMV::rank,
                "KokkosBlas::scal: "
                "R and X must have the same rank.");
  static_assert(RMV::rank == 1 || RMV::rank == 2,
                "KokkosBlas::scal: "
                "RMV and XMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if (X.extent(0) != R.extent(0) || X.extent(1) != R.extent(1)) {
    std::ostringstream os;
    os << "KokkosBlas::scal: Dimensions of R and X do not match: "
       << "R: " << R.extent(0) << " x " << R.extent(1) << ", X: " << X.extent(0) << " x " << X.extent(1);
    KokkosKernels::Impl::throw_runtime_exception(os.str());
  }

  using UnifiedRLayout = typename KokkosKernels::Impl::GetUnifiedLayout<RMV>::array_layout;
  using UnifiedXLayout = typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<XMV, UnifiedRLayout>::array_layout;

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.  AV may be either a rank-1 View, or a scalar
  // value.
  using RMV_Internal = Kokkos::View<typename RMV::non_const_data_type, UnifiedRLayout, typename RMV::device_type,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using XMV_Internal = Kokkos::View<typename XMV::const_data_type, UnifiedXLayout, typename XMV::device_type,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using AV_Internal  = typename KokkosKernels::Impl::GetUnifiedScalarViewType<AV, XMV_Internal, true>::type;

  RMV_Internal R_internal = R;
  AV_Internal a_internal  = a;
  XMV_Internal X_internal = X;

  Impl::Scal<execution_space, RMV_Internal, AV_Internal, XMV_Internal>::scal(space, R_internal, a_internal, X_internal);
}

/// \brief Computes R := alpha*X
///
/// This function is non-blocking and thread-safe
/// The kernel is executed in the default stream/queue
/// associated with the execution space of YMV.
///
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization. It must have
///   the same rank as RMV.
/// \tparam AV 1-D or 2-D Kokkos::View specialization.
///
/// \param R [in/out] view of type RMV in which the results will be stored.
/// \param a [in] view of type AV, scaling parameter for X.
/// \param X [in] input view of type XMV.
template <class RMV, class AV, class XMV>
void scal(const RMV& R, const AV& a, const XMV& X) {
  scal(typename RMV::execution_space{}, R, a, X);
}

///
/// Serial Scale
///

struct SerialScale {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const ScalarType alpha, const AViewType& A) {
    return Impl::SerialScaleInternal::invoke(A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(), A.stride_1());
  }
};

///
/// Team Scale
///

template <typename MemberType>
struct TeamScale {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A) {
    return Impl::TeamScaleInternal::invoke(member, A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(),
                                           A.stride_1());
  }
};

///
/// TeamVector Scale
///

template <typename MemberType>
struct TeamVectorScale {
  template <typename ScalarType, typename AViewType>
  KOKKOS_INLINE_FUNCTION static int invoke(const MemberType& member, const ScalarType alpha, const AViewType& A) {
    return Impl::TeamVectorScaleInternal::invoke(member, A.extent(0), A.extent(1), alpha, A.data(), A.stride_0(),
                                                 A.stride_1());
  }
};

}  // namespace KokkosBlas

#endif
