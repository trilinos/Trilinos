// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_SCAL_UNIFICATION_HPP_
#define KOKKOSBLAS1_SCAL_UNIFICATION_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosKernels_helpers.hpp>

namespace KokkosBlas {
namespace Impl {

// Utilities for unifying a scal coefficient (alpha) into a canonical type
// accepted by the Scal spec/functor layer.
//
// Valid input types for alpha:
//   - scalar                      → scalar (pass through)
//   - rank-0 host-accessible View → extract value on host, use as scalar
//   - rank-0 device-only View     → rank-0 unmanaged View (read in kernel)
//   - rank-1 View (2-D X only)    → rank-1 unmanaged View (one alpha per col)
//
// PreferredScalar should be the value_type of the vector being scaled.
// PreferredLayout should be the layout of the X vector for ETI-friendly dispatch.

// General case: Coeff is not a View → use scalar
// ExecSpace is the kernel execution space; it is used as the device type for unified View
// coefficients so the resulting type matches the ETI declarations (which use EXEC_SPACE directly).
template <typename ExecSpace, typename Coeff, typename PreferredScalar, typename PreferredLayout,
          bool isView = Kokkos::is_view_v<Coeff>>
struct UnifiedScalCoeff {
  using type = PreferredScalar;
};

// Specialization: Coeff is a View
template <typename ExecSpace, typename Coeff, typename PreferredScalar, typename PreferredLayout>
struct UnifiedScalCoeff<ExecSpace, Coeff, PreferredScalar, PreferredLayout, true> {
  static constexpr bool IsRank0 = (int)Coeff::rank == 0;
  static constexpr bool IsHostAccessible =
      Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename Coeff::memory_space>::accessible;
  static constexpr bool IsExecAccessible =
      Kokkos::SpaceAccessibility<ExecSpace, typename Coeff::memory_space>::accessible;

  // Preserves rank: const_data_type of rank-0 is `const scalar`,
  //                 const_data_type of rank-1 is `const scalar*`.
  // Use ExecSpace (not Coeff::device_type) so the resulting type matches the ETI
  // declarations, which are parameterised on the execution space directly.
  using UnifiedViewType =
      Kokkos::View<typename Coeff::const_data_type,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<Coeff, PreferredLayout>::array_layout,
                   ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  // rank-0 + host-only accessibility → scalar; otherwise → rank-preserving unmanaged View
  using type = std::conditional_t<IsRank0 && IsHostAccessible && !IsExecAccessible, PreferredScalar, UnifiedViewType>;
};

template <typename UnifiedCoeff, typename Coeff>
UnifiedCoeff unifyScalCoeff(const Coeff& coeff) {
  if constexpr (Kokkos::is_view_v<Coeff> && !Kokkos::is_view_v<UnifiedCoeff>) {
    // rank-0 host-accessible View → extract scalar value on host
    static_assert((int)Coeff::rank == 0, "unifyScalCoeff: Coeff must be rank-0 in this case");
    static_assert(
        Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename Coeff::memory_space>::accessible,
        "unifyScalCoeff: Coeff must be host-accessible in this case");
    return coeff();
  } else if constexpr (Kokkos::is_view_v<Coeff>) {
    // rank-0 device View or rank-1 View → unmanaged shallow copy
    return KokkosKernels::Impl::unificationCast<UnifiedCoeff>(coeff);
  } else {
    // scalar
    return coeff;
  }
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_SCAL_UNIFICATION_HPP_
