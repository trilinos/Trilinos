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

#ifndef KOKKOSBATCHED_TBSV_SERIAL_INTERNAL_HPP_
#define KOKKOSBATCHED_TBSV_SERIAL_INTERNAL_HPP_

/// \author Yuuichi Asahi (yuuichi.asahi@cea.fr)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {

///
/// Serial Internal Impl
/// ====================

///
/// Lower, Non-Transpose
///

template <typename AlgoType>
struct SerialTbsvInternalLower {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int an, const ValueType *KOKKOS_RESTRICT A,
                                           const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTbsvInternalLower<Algo::Tbsv::Unblocked>::invoke(
    const bool use_unit_diag, const int an, const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = 0; j < an; ++j) {
    if (x[j * xs0] != static_cast<ValueType>(0)) {
      if (!use_unit_diag) x[j * xs0] = x[j * xs0] / A[0 + j * as1];

      auto temp = x[j * xs0];
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = j + 1; i < Kokkos::min(an, j + k + 1); ++i) {
        x[i * xs0] = x[i * xs0] - temp * A[(i - j) * as0 + j * as1];
      }
    }
  }

  return 0;
}

///
/// Lower, Transpose
///

template <typename AlgoType>
struct SerialTbsvInternalLowerTranspose {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const bool do_conj, const int an,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTbsvInternalLowerTranspose<Algo::Tbsv::Unblocked>::invoke(
    const bool use_unit_diag, const bool do_conj, const int an, const ValueType *KOKKOS_RESTRICT A, const int as0,
    const int as1,
    /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = an - 1; j >= 0; --j) {
    auto temp = x[j * xs0];

    if (do_conj) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = Kokkos::min(an - 1, j + k); i > j; --i) {
        temp -= Kokkos::ArithTraits<ValueType>::conj(A[(i - j) * as0 + j * as1]) * x[i * xs0];
      }
      if (!use_unit_diag) temp = temp / Kokkos::ArithTraits<ValueType>::conj(A[0 + j * as1]);
    } else {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = Kokkos::min(an - 1, j + k); i > j; --i) {
        temp -= A[(i - j) * as0 + j * as1] * x[i * xs0];
      }
      if (!use_unit_diag) temp = temp / A[0 + j * as1];
    }
    x[j * xs0] = temp;
  }

  return 0;
}

///
/// Upper, Non-Transpose
///

template <typename AlgoType>
struct SerialTbsvInternalUpper {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const int an, const ValueType *KOKKOS_RESTRICT A,
                                           const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTbsvInternalUpper<Algo::Tbsv::Unblocked>::invoke(
    const bool use_unit_diag, const int an, const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
    /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = an - 1; j >= 0; --j) {
    if (x[j * xs0] != 0) {
      if (!use_unit_diag) x[j * xs0] = x[j * xs0] / A[k * as0 + j * as1];

      auto temp = x[j * xs0];
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = j - 1; i >= Kokkos::max(0, j - k); --i) {
        x[i * xs0] = x[i * xs0] - temp * A[(k - j + i) * as0 + j * as1];
      }
    }
  }

  return 0;
}

///
/// Upper, Transpose
///

template <typename AlgoType>
struct SerialTbsvInternalUpperTranspose {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION static int invoke(const bool use_unit_diag, const bool do_conj, const int an,
                                           const ValueType *KOKKOS_RESTRICT A, const int as0, const int as1,
                                           /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k);
};

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int SerialTbsvInternalUpperTranspose<Algo::Tbsv::Unblocked>::invoke(
    const bool use_unit_diag, const bool do_conj, const int an, const ValueType *KOKKOS_RESTRICT A, const int as0,
    const int as1,
    /**/ ValueType *KOKKOS_RESTRICT x, const int xs0, const int k) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = 0; j < an; j++) {
    auto temp = x[j * xs0];
    if (do_conj) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = Kokkos::max(0, j - k); i < j; ++i) {
        temp -= Kokkos::ArithTraits<ValueType>::conj(A[(i + k - j) * as0 + j * as1]) * x[i * xs0];
      }
      if (!use_unit_diag) temp = temp / Kokkos::ArithTraits<ValueType>::conj(A[k * as0 + j * as1]);
    } else {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      for (int i = Kokkos::max(0, j - k); i < j; ++i) {
        temp -= A[(i + k - j) * as0 + j * as1] * x[i * xs0];
      }
      if (!use_unit_diag) temp = temp / A[k * as0 + j * as1];
    }
    x[j * xs0] = temp;
  }

  return 0;
}

}  // namespace KokkosBatched

#endif  // KOKKOSBATCHED_TBSV_SERIAL_INTERNAL_HPP_
