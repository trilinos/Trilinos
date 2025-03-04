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

#ifndef KOKKOSBATCHED_LACGV_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_LACGV_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ====================

struct SerialLacgvInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int xn, ValueType *KOKKOS_RESTRICT x, const int xs0);

  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int xn, Kokkos::complex<ValueType> *KOKKOS_RESTRICT x, const int xs0);
};

// Real specialization (no op)
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialLacgvInternal::invoke(const int /*xn*/, ValueType *KOKKOS_RESTRICT /*x*/,
                                                       const int /*xs0*/) {
  return 0;
}

// Complex specialization
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialLacgvInternal::invoke(const int xn, Kokkos::complex<ValueType> *KOKKOS_RESTRICT x,
                                                       const int xs0) {
  for (int i = 0; i < xn; i++) {
    x[i * xs0] = Kokkos::conj(x[i * xs0]);
  }
  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_LACGV_SERIAL_INTERNAL_HPP_
