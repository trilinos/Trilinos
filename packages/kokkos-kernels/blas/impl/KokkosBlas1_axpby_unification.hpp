// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBLAS1_AXPBY_UNIFICATION_HPP_
#define KOKKOSBLAS1_AXPBY_UNIFICATION_HPP_

#include <Kokkos_Core.hpp>
#include <KokkosKernels_helpers.hpp>

namespace KokkosBlas {
namespace Impl {

template <class T>
  requires(Kokkos::is_view_v<T>)
constexpr KOKKOS_INLINE_FUNCTION bool isRank1View() {
  return (int)T::rank == 1;
}

template <class T>
  requires(!Kokkos::is_view_v<T>)
constexpr KOKKOS_INLINE_FUNCTION bool isRank1View() {
  return false;
}

// Utilities for unifying an axpby coefficient (a or b) in rank-1 X/Y case.
// When Coeff is a rank-0 HostSpace View, read its value to get a scalar.
// Then use the preferred scalar type (if coeff is a scalar) or
// rank-1 view with a given layout (if coeff is a View).
//
// PreferredScalar should be the value type of the vector that this coefficient will be applied to.
// PreferredLayout should be the layout of the X vector, since in the best case for ETI, A, X, B, Y will
// all have the same layout.

// General case for isView == false
template <typename Coeff, typename PreferredScalar, typename PreferredLayout, bool isView = Kokkos::is_view_v<Coeff>>
struct UnifiedAxpbyCoeff {
  using type = PreferredScalar;
};

// Specialization for isView == true
template <typename Coeff, typename PreferredScalar, typename PreferredLayout>
struct UnifiedAxpbyCoeff<Coeff, PreferredScalar, PreferredLayout, true> {
  static constexpr bool IsRank0 = (int)Coeff::rank == 0;
  static constexpr bool IsHostAccessible =
      Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename Coeff::memory_space>::accessible;

  using UnifiedViewType =
      Kokkos::View<typename Coeff::const_value_type*,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<Coeff, PreferredLayout>::array_layout,
                   typename Coeff::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using type = std::conditional_t<IsRank0 && IsHostAccessible, PreferredScalar, UnifiedViewType>;
};

template <typename UnifiedCoeff, typename Coeff>
UnifiedCoeff unifyAxpbyCoeff(const Coeff& coeff) {
  if constexpr (Kokkos::is_view_v<Coeff> && !Kokkos::is_view_v<UnifiedCoeff>) {
    // Access the scalar value of rank-0, HostSpace coeff
    static_assert((int)Coeff::rank == 0, "unifyAxpbyCoeff: in this case the rank of Coeff must be 0");
    static_assert(
        Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename Coeff::memory_space>::accessible,
        "unifyAxpbyCoeff: in this case Coeff needs to be a host accessible View");
    return coeff();
  } else if constexpr (isRank1View<Coeff>()) {
    // Directly convert to unified type
    return KokkosKernels::Impl::unificationCast<UnifiedCoeff>(coeff);
  } else if constexpr (isRank1View<UnifiedCoeff>()) {
    // UnifiedCoeff is a rank-1 View but Coeff is rank-0.
    // Convert rank-0 Coeff to rank-1 UnifiedCoeff
    return UnifiedCoeff(coeff.data(), 1);
  } else {
    // Coeff is a scalar
    return coeff;
  }
}

// For AxpbyMv, the unified coefficient type is the same as above except that
// we do not change rank-0 Views to rank-1.
template <typename Coeff, typename PreferredScalar, typename PreferredLayout, bool isView = Kokkos::is_view_v<Coeff>>
struct UnifiedAxpbyMvCoeff {
  using type = PreferredScalar;
};

// Specialization for isView == true
template <typename Coeff, typename PreferredScalar, typename PreferredLayout>
struct UnifiedAxpbyMvCoeff<Coeff, PreferredScalar, PreferredLayout, true> {
  static constexpr bool IsRank0 = (int)Coeff::rank == 0;
  static constexpr bool IsHostAccessible =
      Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename Coeff::memory_space>::accessible;

  using UnifiedViewType =
      Kokkos::View<typename Coeff::const_data_type,
                   typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<Coeff, PreferredLayout>::array_layout,
                   typename Coeff::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  using type = std::conditional_t<IsRank0 && IsHostAccessible, PreferredScalar, UnifiedViewType>;
};

template <typename UnifiedCoeff, typename Coeff>
UnifiedCoeff unifyAxpbyMvCoeff(const Coeff& coeff) {
  if constexpr (Kokkos::is_view_v<Coeff> && !Kokkos::is_view_v<UnifiedCoeff>) {
    // Access the scalar value of rank-0, HostSpace coeff
    static_assert((int)Coeff::rank == 0, "unifyAxpbyCoeff: in this case the rank of Coeff must be 0");
    static_assert(
        Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, typename Coeff::memory_space>::accessible,
        "unifyAxpbyCoeff: in this case Coeff needs to be a host accessible View");
    return coeff();
  } else {
    return KokkosKernels::Impl::unificationCast<UnifiedCoeff>(coeff);
  }
}

template <typename Coeff, bool isView = Kokkos::is_view_v<Coeff>>
struct CoeffScalarType {
  using type = Coeff;
};

template <typename Coeff>
struct CoeffScalarType<Coeff, true> {
  using type = typename Coeff::non_const_value_type;
};

// Extract the value from c from any of these cases:
// - c is a scalar
// - c is a rank-0 device view
// - c is a rank-1 device view (with one coefficient per column)
template <typename Coeff>
KOKKOS_INLINE_FUNCTION typename CoeffScalarType<Coeff>::type getCoefficient(const Coeff& c, int column = 0) {
  if constexpr (Kokkos::is_view_v<Coeff>) {
    if constexpr (isRank1View<Coeff>()) {
      return c(column);
    } else {
      // rank-0
      return c();
    }
  } else {
    // scalar
    return c;
  }
}

// Version of getCoefficient that supports c being a rank-1, extent 1 View (where column can be anything).
// This is the "slow" version because it uses a runtime branch
template <typename Coeff>
KOKKOS_INLINE_FUNCTION typename CoeffScalarType<Coeff>::type getCoefficientSlow(const Coeff& c, int column = 0) {
  if constexpr (Kokkos::is_view_v<Coeff>) {
    if constexpr (isRank1View<Coeff>()) {
      if (c.extent_int(0) == 1)
        return c(0);
      else
        return c(column);
    } else {
      // rank-0
      return c();
    }
  } else {
    // scalar
    return c;
  }
}

// Get the number of axpby coefficients represented by c in any of these cases:
// - c is a scalar
// - c is a rank-0 view
// - c is a rank-1 view
template <typename Coeff>
KOKKOS_INLINE_FUNCTION int getNumCoefficients(const Coeff& c) {
  if constexpr (Kokkos::is_view_v<Coeff>) {
    if constexpr (isRank1View<Coeff>()) {
      return c.extent(0);
    } else {
      // rank-0
      return 1;
    }
  } else {
    // scalar
    return 1;
  }
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS1_AXPBY_UNIFICATION_HPP_
