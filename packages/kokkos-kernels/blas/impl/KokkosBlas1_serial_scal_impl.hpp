// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS1_SERIAL_SCAL_IMPL_HPP_
#define KOKKOSBLAS1_SERIAL_SCAL_IMPL_HPP_

#include <Kokkos_Core.hpp>

namespace KokkosBlas {
namespace Impl {

///
/// Serial Internal Impl
/// ====================
struct SerialScaleInternal {
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < m; ++i) A[i * as0] *= alpha;

    return 0;
  }

  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int m, const int n, const ScalarType alpha,
                                           /* */ ValueType *KOKKOS_RESTRICT A, const int as0, const int as1) {
    if (as0 > as1)
      for (int i = 0; i < m; ++i) invoke(n, alpha, A + i * as0, as1);
    else
      for (int j = 0; j < n; ++j) invoke(m, alpha, A + j * as1, as0);

    return 0;
  }
};

}  // namespace Impl
}  // namespace KokkosBlas

#endif
