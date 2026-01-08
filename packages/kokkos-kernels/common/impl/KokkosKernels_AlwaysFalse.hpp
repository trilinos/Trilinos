// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSKERNELS_ALWAYSFALSE_HPP
#define KOKKOSKERNELS_ALWAYSFALSE_HPP

namespace KokkosKernels::Impl {

// for use in static asserts
template <typename...>
inline constexpr bool always_false_v = false;

}  // namespace KokkosKernels::Impl

#endif  // KOKKOSKERNELS_ALWAYSFALSE_HPP
