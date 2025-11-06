// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSBLAS_SERIAL_NRM2_HPP_
#define KOKKOSBLAS_SERIAL_NRM2_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

///
/// Serial Internal Impl
/// ====================
template <typename ValueType>
KOKKOS_INLINE_FUNCTION static typename Kokkos::Details::InnerProductSpaceTraits<ValueType>::mag_type serial_nrm2(
    const int m, const ValueType *KOKKOS_RESTRICT X, const int xs0) {
  using IPT       = Kokkos::Details::InnerProductSpaceTraits<ValueType>;
  using norm_type = typename IPT::mag_type;

  norm_type nrm = KokkosKernels::ArithTraits<norm_type>::zero();

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int i = 0; i < m; ++i) nrm += IPT::norm(IPT::dot(X[i * xs0], X[i * xs0]));

  return KokkosKernels::ArithTraits<norm_type>::sqrt(nrm);
}

template <typename ValueType>
KOKKOS_INLINE_FUNCTION static void serial_nrm2(
    const int m, const int n, const ValueType *KOKKOS_RESTRICT X, const int xs0, const int xs1,
    typename Kokkos::Details::InnerProductSpaceTraits<ValueType>::mag_type *KOKKOS_RESTRICT R, const int ys0) {
  for (int vecIdx = 0; vecIdx < n; ++vecIdx) R[vecIdx * ys0] = serial_nrm2(m, X + vecIdx * xs1, xs0);

  return;
}

}  // namespace Impl
}  // namespace KokkosBlas

#endif  // KOKKOSBLAS_SERIAL_NRM2_HPP_
