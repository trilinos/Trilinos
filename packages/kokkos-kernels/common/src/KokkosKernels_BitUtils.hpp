// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_BITUTILS_HPP
#define KOKKOSKERNELS_BITUTILS_HPP
#include "Kokkos_Core.hpp"

#include <type_traits>

namespace KokkosKernels::Impl {

template <class T>
KOKKOS_FUNCTION std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> pop_count(T x) {
  return Kokkos::Experimental::popcount_builtin<std::make_unsigned_t<T>>(x);
}

template <class T>
KOKKOS_FUNCTION std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<T, bool>, int> least_set_bit(T x) {
  return Kokkos::Experimental::countr_zero_builtin<std::make_unsigned_t<T>>(x) + 1;
}

}  // namespace KokkosKernels::Impl

#endif
