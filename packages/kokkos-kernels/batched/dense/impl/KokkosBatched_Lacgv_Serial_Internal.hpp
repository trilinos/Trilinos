// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
