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
      typename T::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      type;
};

template <class T, class TX>
struct GetUnifiedScalarViewType<T, TX, true, true> {
  typedef Kokkos::View<
      typename T::const_value_type*,
      typename KokkosKernels::Impl::GetUnifiedLayoutPreferring<T, typename TX::array_layout>::array_layout,
      typename T::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
      type;
};

template <class... Ts>
struct are_integral : std::bool_constant<((std::is_integral_v<Ts> || std::is_enum_v<Ts>)&&...)> {};

template <class... Ts>
inline constexpr bool are_integral_v = are_integral<Ts...>::value;

}  // namespace Impl
}  // namespace KokkosKernels
#endif
