// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_PTTRF_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_PTTRF_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {
namespace Impl {
template <typename AlgoType>
struct SerialPttrfInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, ValueType *KOKKOS_RESTRICT d, const int ds0,
                                           ValueType *KOKKOS_RESTRICT e, const int es0);

  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, ValueType *KOKKOS_RESTRICT d, const int ds0,
                                           Kokkos::complex<ValueType> *KOKKOS_RESTRICT e, const int es0);
};

///
/// Real matrix
///

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPttrfInternal<Algo::Pttrf::Unblocked>::invoke(
    const int n, ValueType *KOKKOS_RESTRICT d, const int ds0, ValueType *KOKKOS_RESTRICT e, const int es0) {
  auto update = [&](const int i) {
    auto ei_tmp = e[i * es0];
    e[i * es0]  = ei_tmp / d[i * ds0];
    d[(i + 1) * ds0] -= e[i * es0] * ei_tmp;
  };

#ifndef NDEBUG
  int info = 0;

  auto check_positive_definitiveness = [&](const int i) {
    return (d[i] <= KokkosKernels::ArithTraits<ValueType>::zero()) ? (i + 1) : 0;
  };
#endif

  // Compute the L*D*L' (or U'*D*U) factorization of A.
  const int i4 = (n - 1) % 4;
  for (int i = 0; i < i4; i++) {
#ifndef NDEBUG
    info = check_positive_definitiveness(i);
    if (info) return info;
#endif

    update(i);
  }  // for (int i = 0; i < i4; i++)

  for (int i = i4; i < n - 4; i += 4) {
#ifndef NDEBUG
    info = check_positive_definitiveness(i);
    if (info) return info;
#endif

    update(i);

#ifndef NDEBUG
    info = check_positive_definitiveness(i + 1);
    if (info) return info;
#endif

    update(i + 1);

#ifndef NDEBUG
    info = check_positive_definitiveness(i + 2);
    if (info) return info;
#endif

    update(i + 2);

#ifndef NDEBUG
    info = check_positive_definitiveness(i + 3);
    if (info) return info;
#endif

    update(i + 3);
  }  // for (int i = i4; i < n-4; 4)

#ifndef NDEBUG
  info = check_positive_definitiveness(n - 1);
  if (info) return info;
#endif

  return 0;
}

///
/// Complex matrix
///

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPttrfInternal<Algo::Pttrf::Unblocked>::invoke(
    const int n, ValueType *KOKKOS_RESTRICT d, const int ds0, Kokkos::complex<ValueType> *KOKKOS_RESTRICT e,
    const int es0) {
  auto update = [&](const int i) {
    auto eir_tmp     = e[i * es0].real();
    auto eii_tmp     = e[i * es0].imag();
    auto f_tmp       = eir_tmp / d[i * ds0];
    auto g_tmp       = eii_tmp / d[i * ds0];
    e[i * es0]       = Kokkos::complex<ValueType>(f_tmp, g_tmp);
    d[(i + 1) * ds0] = d[(i + 1) * ds0] - f_tmp * eir_tmp - g_tmp * eii_tmp;
  };

#ifndef NDEBUG
  int info = 0;

  auto check_positive_definitiveness = [&](const int i) {
    return (d[i] <= KokkosKernels::ArithTraits<ValueType>::zero()) ? (i + 1) : 0;
  };
#endif

  // Compute the L*D*L' (or U'*D*U) factorization of A.
  const int i4 = (n - 1) % 4;
  for (int i = 0; i < i4; i++) {
#ifndef NDEBUG
    info = check_positive_definitiveness(i);
    if (info) return info;
#endif

    update(i);
  }  // for (int i = 0; i < i4; i++)

  for (int i = i4; i < n - 4; i += 4) {
#ifndef NDEBUG
    info = check_positive_definitiveness(i);
    if (info) return info;
#endif

    update(i);

#ifndef NDEBUG
    info = check_positive_definitiveness(i + 1);
    if (info) return info;
#endif

    update(i + 1);

#ifndef NDEBUG
    info = check_positive_definitiveness(i + 2);
    if (info) return info;
#endif

    update(i + 2);

#ifndef NDEBUG
    info = check_positive_definitiveness(i + 3);
    if (info) return info;
#endif

    update(i + 3);
  }  // for (int i = i4; i < n-4; 4)

#ifndef NDEBUG
  info = check_positive_definitiveness(n - 1);
  if (info) return info;
#endif

  return 0;
}

}  // namespace Impl
}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PTTRF_SERIAL_INTERNAL_HPP_
