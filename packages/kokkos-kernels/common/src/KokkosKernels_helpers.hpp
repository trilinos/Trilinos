// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSKERNELS_HELPERS_HPP_
#define KOKKOSKERNELS_HELPERS_HPP_

#include "KokkosKernels_config.h"           // KOKKOSKERNELS_INST_LAYOUTLEFT, KOKKOSKERNELS_INST_LAYOUTRIGHT
#include "KokkosKernels_default_types.hpp"  // default_layout

#include <type_traits>

namespace KokkosKernels {
namespace Impl {

// Unify Layout of a View to PreferredLayoutType if possible
// (either matches already, or is rank-0/rank-1 and contiguous)
// Used to reduce number of code instantiations.
template <class ViewType, class PreferredLayoutType>
struct GetUnifiedLayoutPreferring {
  using array_layout =
      typename std::conditional<((ViewType::rank == 1) &&
                                 !std::is_same_v<typename ViewType::array_layout, Kokkos::LayoutStride>) ||
                                    (ViewType::rank == 0),
                                PreferredLayoutType, typename ViewType::array_layout>::type;
};

template <class ViewType>
struct GetUnifiedLayout {
  using array_layout = typename GetUnifiedLayoutPreferring<ViewType, default_layout>::array_layout;
};

template <class T, class TX, bool do_const, bool isView = Kokkos::is_view<T>::value>
struct GetUnifiedScalarViewType {
  typedef typename TX::non_const_value_type type;
};

template <class T, class TX>
struct GetUnifiedScalarViewType<T, TX, false, true> {
  typedef Kokkos::View<
      typename T::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<T, typename TX::array_layout>::array_layout,
      typename T::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      type;
};

template <class T, class TX>
struct GetUnifiedScalarViewType<T, TX, true, true> {
  typedef Kokkos::View<
      typename T::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<T, typename TX::array_layout>::array_layout,
      typename T::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      type;
};

// GetUnifiedScalarExecSpaceViewType: like the above, unifies a parameter type T that could be either a scalar or View.
// But if it's a View, we use ExecSpace as the unified device type (instead of T::device_type as with
// GetUnifiedScalarViewType above)
template <class ExecSpace, class T, class TX, bool do_const, bool isView = Kokkos::is_view<T>::value>
struct GetUnifiedScalarExecSpaceViewType {
  typedef typename TX::non_const_value_type type;
};

template <class ExecSpace, class T, class TX>
struct GetUnifiedScalarExecSpaceViewType<ExecSpace, T, TX, false, true> {
  typedef Kokkos::View<
      typename T::non_const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<T, typename TX::array_layout>::array_layout, ExecSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      type;
};

template <class ExecSpace, class T, class TX>
struct GetUnifiedScalarExecSpaceViewType<ExecSpace, T, TX, true, true> {
  typedef Kokkos::View<
      typename T::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<T, typename TX::array_layout>::array_layout, ExecSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      type;
};

// Utility for shallow-copying a view v to a different type, where
// direct assignment may or may not be allowed by Kokkos.
//
// Some situations where direct assignment is not allowed (and this helper is needed):
// - V1 is HostSpace, and V2 is any shared space.
// - V1 is Serial and V2 is Device<Serial, CudaUVMSpace>
//
// In some kernels, unified parameters can be scalars so this also works with V1 and V2 being convertible scalar types.
template <typename V1, typename V2>
V1 unificationCast(const V2& v) {
  if constexpr (!Kokkos::is_view_v<V1>) {
    return v;
  } else {
    if constexpr (std::is_same_v<typename V1::array_layout, typename V2::array_layout>) {
      return V1(v.data(), v.layout());
    }
    // If layout types aren't the same, first convert to an intermediate type with the target layout
    Kokkos::View<typename V2::data_type, typename V1::array_layout, typename V2::device_type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        temp = v;
    // Then we can we use View assignment to change to V1's device type and traits
    return V1(temp.data(), temp.layout());
  }
}

template <class... Ts>
struct are_integral : std::bool_constant<((std::is_integral_v<Ts> || std::is_enum_v<Ts>)&&...)> {};

template <class... Ts>
inline constexpr bool are_integral_v = are_integral<Ts...>::value;

}  // namespace Impl
}  // namespace KokkosKernels
#endif
