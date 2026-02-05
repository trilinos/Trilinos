// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBATCHED_IAMAX_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_IAMAX_SERIAL_INTERNAL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include <Kokkos_Core.hpp>
#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
namespace Impl {

///
/// Serial Internal Impl
/// ========================

struct SerialIamaxInternal {
  template <typename IndexType, typename ValueType>
  KOKKOS_INLINE_FUNCTION static IndexType invoke(const int n, const ValueType *KOKKOS_RESTRICT x, const int xs0);
};

template <typename IndexType, typename ValueType>
KOKKOS_INLINE_FUNCTION IndexType SerialIamaxInternal::invoke(const int n, const ValueType *KOKKOS_RESTRICT x,
                                                             const int xs0) {
  using ats      = typename KokkosKernels::ArithTraits<ValueType>;
  using RealType = typename ats::mag_type;

  RealType amax  = Kokkos::abs(x[0 * xs0]);
  IndexType imax = 0;

  for (IndexType i = 1; i < static_cast<IndexType>(n); ++i) {
    const RealType abs_x_i = Kokkos::abs(x[i * xs0]);
    if (abs_x_i > amax) {
      amax = abs_x_i;
      imax = i;
    }
  }

  return imax;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_IAMAX_SERIAL_INTERNAL_HPP_
