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
#ifndef KOKKOSBATCHED_PTTRS_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_PTTRS_SERIAL_INTERNAL_HPP_

#include <KokkosBatched_Util.hpp>

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

namespace KokkosBatched {

template <typename ArgUplo, typename AlgoType>
struct SerialPttrsInternal {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, const ValueType *KOKKOS_RESTRICT d, const int ds0,
                                           const ValueType *KOKKOS_RESTRICT e, const int es0,
                                           ValueType *KOKKOS_RESTRICT b, const int bs0);

  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const int n, const ValueType *KOKKOS_RESTRICT d, const int ds0,
                                           const Kokkos::complex<ValueType> *KOKKOS_RESTRICT e, const int es0,
                                           Kokkos::complex<ValueType> *KOKKOS_RESTRICT b, const int bs0);
};

///
/// Real matrix
///

template <typename ArgUplo, typename AlgoType>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPttrsInternal<ArgUplo, AlgoType>::invoke(
    const int n, const ValueType *KOKKOS_RESTRICT d, const int ds0, const ValueType *KOKKOS_RESTRICT e, const int es0,
    ValueType *KOKKOS_RESTRICT b, const int bs0) {
  // Solve A * X = B using the factorization L * D * L**T
  for (int i = 1; i < n; i++) {
    b[i * bs0] -= e[(i - 1) * es0] * b[(i - 1) * bs0];
  }

  b[(n - 1) * bs0] /= d[(n - 1) * ds0];

  for (int i = n - 2; i >= 0; i--) {
    b[i * bs0] = b[i * bs0] / d[i * ds0] - b[(i + 1) * bs0] * e[i * es0];
  }

  return 0;
}

///
/// Complex matrix
///

template <typename ArgUplo, typename AlgoType>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialPttrsInternal<ArgUplo, AlgoType>::invoke(
    const int n, const ValueType *KOKKOS_RESTRICT d, const int ds0, const Kokkos::complex<ValueType> *KOKKOS_RESTRICT e,
    const int es0, Kokkos::complex<ValueType> *KOKKOS_RESTRICT b, const int bs0) {
  // Solve A * X = B using the factorization L * D * L**H
  for (int i = 1; i < n; i++) {
    auto tmp_e = std::is_same_v<ArgUplo, Uplo::Upper> ? Kokkos::conj(e[(i - 1) * es0]) : e[(i - 1) * es0];
    b[i * bs0] -= tmp_e * b[(i - 1) * bs0];
  }

  b[(n - 1) * bs0] /= d[(n - 1) * ds0];

  for (int i = n - 2; i >= 0; i--) {
    auto tmp_e = std::is_same_v<ArgUplo, Uplo::Lower> ? Kokkos::conj(e[i * es0]) : e[i * es0];
    b[i * bs0] = b[i * bs0] / d[i * ds0] - b[(i + 1) * bs0] * tmp_e;
  }

  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_PTTRS_SERIAL_INTERNAL_HPP_
