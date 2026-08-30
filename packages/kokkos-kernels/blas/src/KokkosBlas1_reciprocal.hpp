// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_RECIPROCAL_HPP_
#define KOKKOSBLAS1_RECIPROCAL_HPP_

#include <KokkosBlas1_reciprocal_spec.hpp>
#include <KokkosKernels_helpers.hpp>
#include <KokkosKernels_Error.hpp>

namespace KokkosBlas {

/// \brief R(i,j) = reciprocal(X(i,j))
///
/// Replace each entry in R with the absolute value (magnitude), of the
/// reciprocal of the corresponding entry in X.
/// This function is non-blocking and thread-safe
///
/// \tparam execution_space a Kokkos execution space
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
///
/// \param space [in] an instance of execution space where the kernel will run
/// \param R [out] a view of type RMV that contains the inverse of the values in
/// X.
/// \param X [in] a view of type XMV that contains the values to invert.
template <class execution_space, class RMV, class XMV>
void reciprocal(const execution_space& space, const RMV& R, const XMV& X) {
  static_assert(Kokkos::is_execution_space_v<execution_space>,
                "KokkosBlas::reciprocal: execution_space must be a valid "
                "Kokkos execition space.");
  static_assert(Kokkos::is_view<RMV>::value,
                "KokkosBlas::reciprocal: "
                "R is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename RMV::memory_space>::accessible,
                "KokkosBlas::reciprocal: RMV must be accessible from execution_space");
  static_assert(Kokkos::is_view<XMV>::value,
                "KokkosBlas::reciprocal: "
                "X is not a Kokkos::View.");
  static_assert(Kokkos::SpaceAccessibility<execution_space, typename XMV::memory_space>::accessible,
                "KokkosBlas::reciprocal: XMV must be accessible from execution_space");
  static_assert(std::is_same<typename RMV::value_type, typename RMV::non_const_value_type>::value,
                "KokkosBlas::reciprocal: R is const.  "
                "It must be nonconst, because it is an output argument "
                "(we have to be able to write to its entries).");
  static_assert(int(RMV::rank) == int(XMV::rank),
                "KokkosBlas::reciprocal: "
                "R and X must have the same rank.");
  static_assert(RMV::rank == 1 || RMV::rank == 2,
                "KokkosBlas::reciprocal: "
                "RMV and XMV must either have rank 1 or rank 2.");

  // Check compatibility of dimensions at run time.
  if constexpr (RMV::rank == 1) {
    if (X.extent(0) != R.extent(0)) {
      std::ostringstream os;
      os << "KokkosBlas::reciprocal (MV): Dimensions of R and X do not match: "
         << "R: " << R.extent(0) << ", X: " << X.extent(0);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  } else if constexpr (RMV::rank == 2) {
    if (X.extent(0) != R.extent(0) || X.extent(1) != R.extent(1)) {
      std::ostringstream os;
      os << "KokkosBlas::reciprocal (MV): Dimensions of R and X do not match: "
         << "R: " << R.extent(0) << " x " << R.extent(1) << ", X: " << X.extent(0) << " x " << X.extent(1);
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
  }

  // Create unmanaged versions of the input Views.  RMV and XMV may be
  // rank 1 or rank 2.
  using PreferredLayout = typename KokkosKernels::Impl::GetUnifiedLayout<XMV>::array_layout;
  using RMV_Internal =
      Kokkos::View<typename RMV::non_const_data_type,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<RMV, PreferredLayout>::array_layout,
                   execution_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  using XMV_Internal = Kokkos::View<typename XMV::const_data_type, PreferredLayout, execution_space,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  RMV_Internal R_internal = KokkosKernels::Impl::unificationCast<RMV_Internal>(R);
  XMV_Internal X_internal = KokkosKernels::Impl::unificationCast<XMV_Internal>(X);

  Impl::Reciprocal<execution_space, RMV_Internal, XMV_Internal>::reciprocal(space, R_internal, X_internal);
}

/// \brief R(i,j) = reciprocal(X(i,j))
///
/// Replace each entry in R with the absolute value (magnitude), of the
/// reciprocal of the corresponding entry in X.
/// This function is non-blocking and thread-safe
/// The kernel is executed in the default stream/queue
/// associated with the execution space of RMV.
///
/// \tparam RMV 1-D or 2-D Kokkos::View specialization.
/// \tparam XMV 1-D or 2-D Kokkos::View specialization.  It must have
///   the same rank as RMV, and its entries must be assignable to
///   those of RMV.
template <class RMV, class XMV>
void reciprocal(const RMV& R, const XMV& X) {
  reciprocal(typename RMV::execution_space{}, R, X);
}
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_RECIPROCAL_HPP_
