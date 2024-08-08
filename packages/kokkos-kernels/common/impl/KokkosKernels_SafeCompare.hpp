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

#ifndef KOKKOSKERNELS_SAFECOMPARE_HPP
#define KOKKOSKERNELS_SAFECOMPARE_HPP

#include "Kokkos_ArithTraits.hpp"

namespace KokkosKernels {
namespace Impl {

/*! \brief t > u

   When comparing signed and unsigned types of the same size, the signed type
   is converted to unsigned which produces strange behavior like int32_t(-1) >
   uint32_t(1) This function casts its arguments to types that can represent
   the full range of both argument types, before comparing.

    Basically this boils down to:
    1. forbidding any comparisons between signed integers and uint64_t,
    since there's no reliable signed integer type larger than 64 bits.
    2. Using a type large enough to represent both sides of a comparison
   otherwise.

    If T and A are ints, and T xor U is signed, choose a signed type large
    enough to represent all values of both T and U

    This function does not protect you from casting an int to a float where that
   value is not representable.
*/
template <typename T, typename U>
KOKKOS_INLINE_FUNCTION constexpr bool safe_gt(const T &t, const U &u) {
  using KT = Kokkos::ArithTraits<T>;
  using KU = Kokkos::ArithTraits<U>;

  // both are integer, but only one is signed
  if constexpr (KT::is_integer && KU::is_integer && (KT::is_signed != KU::is_signed)) {
    // how wide the signed type would need to be to hold T and U
    constexpr size_t t_width = KT::is_signed ? sizeof(T) : 2 * sizeof(T);
    constexpr size_t u_width = KU::is_signed ? sizeof(U) : 2 * sizeof(U);

    // compare using the max width
    constexpr size_t width = KOKKOSKERNELS_MACRO_MAX(t_width, u_width);
    if constexpr (width == 1) {
      return int8_t(t) > int8_t(u);
    } else if constexpr (width == 2) {
      return int16_t(t) > int16_t(u);
    } else if constexpr (width == 4) {
      return int32_t(t) > int32_t(u);
    } else if constexpr (width == 8) {
      return int64_t(t) > int64_t(u);
    } else {
      static_assert(std::is_same_v<T, U>, "no safe way to compare types");
    }
  } else {
    // use whatever the default comparison rules are
    return t > u;
  }

  // CUDA 11.2 issues a spurious missing return warning
  return false;
}

}  // namespace Impl
}  // namespace KokkosKernels

#endif  // KOKKOSKERNELS_SAFECOMPARE_HPP