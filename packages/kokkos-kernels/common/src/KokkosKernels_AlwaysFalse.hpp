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

#ifndef KOKKOSKERNELS_ALWAYSFALSE_HPP
#define KOKKOSKERNELS_ALWAYSFALSE_HPP

#include <type_traits>

/*! \file KokkosKernels_AlwaysFalse.hpp
    \brief A convenience type to be used in a static_assert that should always
   fail
*/

namespace KokkosKernels {
namespace Impl {

template <typename T>
using always_false = std::false_type;

template <typename T>
inline constexpr bool always_false_v = always_false<T>::value;

}  // namespace Impl
}  // namespace KokkosKernels

#endif  // KOKKOSKERNELS_ALWAYSFALSE_HPP
