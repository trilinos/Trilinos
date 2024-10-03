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
#ifndef __KOKKOSBLAS_INNER_MULTIPLE_DOT_PRODUCT_SERIAL_IMPL_HPP__
#define __KOKKOSBLAS_INNER_MULTIPLE_DOT_PRODUCT_SERIAL_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBlas {
namespace Impl {

struct OpID {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION ValueType operator()(ValueType v) const {
    return v;
  }
};

struct OpConj {
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION ValueType operator()(ValueType v) const {
    using KAT = Kokkos::ArithTraits<ValueType>;
    return KAT::conj(v);
  }
};

template <int mb>
struct InnerMultipleDotProduct {
  const int _as0, _as1, _xs0, _ys0;

  KOKKOS_INLINE_FUNCTION
  InnerMultipleDotProduct(const int as0, const int as1, const int xs0, const int ys0)
      : _as0(as0), _as1(as1), _xs0(xs0), _ys0(ys0) {}

  template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ScalarType alpha, const ValueAType *KOKKOS_RESTRICT A,
                                           const ValueXType *KOKKOS_RESTRICT x, const int n,
                                           ValueYType *KOKKOS_RESTRICT y);

  template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ScalarType alpha, const ValueAType *KOKKOS_RESTRICT A,
                                           const ValueXType *KOKKOS_RESTRICT x, const int m, const int n,
                                           ValueYType *KOKKOS_RESTRICT y);
};

///
/// Dot Product for GEMV
/// ====================

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<5>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int n,
                                                                     ValueYType *KOKKOS_RESTRICT y) {
  if (n <= 0) return 0;

  const int i0 = 0 * _as0, i1 = 1 * _as0, i2 = 2 * _as0, i3 = 3 * _as0, i4 = 4 * _as0;

  // unroll by rows
  ValueYType y_0 = 0, y_1 = 0, y_2 = 0, y_3 = 0, y_4 = 0;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = 0; j < n; ++j) {
    const int jj         = j * _as1;
    const ValueXType x_j = x[j * _xs0];

    y_0 += op(A[i0 + jj]) * x_j;
    y_1 += op(A[i1 + jj]) * x_j;
    y_2 += op(A[i2 + jj]) * x_j;
    y_3 += op(A[i3 + jj]) * x_j;
    y_4 += op(A[i4 + jj]) * x_j;
  }

  y[0 * _ys0] += alpha * y_0;
  y[1 * _ys0] += alpha * y_1;
  y[2 * _ys0] += alpha * y_2;
  y[3 * _ys0] += alpha * y_3;
  y[4 * _ys0] += alpha * y_4;

  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<4>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int n,
                                                                     ValueYType *KOKKOS_RESTRICT y) {
  if (!n) return 0;
  OpA op;

  const int i0 = 0 * _as0, i1 = 1 * _as0, i2 = 2 * _as0, i3 = 3 * _as0;

  // unroll by rows
  ValueYType y_0 = 0, y_1 = 0, y_2 = 0, y_3 = 0;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = 0; j < n; ++j) {
    const int jj         = j * _as1;
    const ValueXType x_j = x[j * _xs0];

    y_0 += op(A[i0 + jj]) * x_j;
    y_1 += op(A[i1 + jj]) * x_j;
    y_2 += op(A[i2 + jj]) * x_j;
    y_3 += op(A[i3 + jj]) * x_j;
  }

  y[0 * _ys0] += alpha * y_0;
  y[1 * _ys0] += alpha * y_1;
  y[2 * _ys0] += alpha * y_2;
  y[3 * _ys0] += alpha * y_3;

  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<3>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int n,
                                                                     ValueYType *KOKKOS_RESTRICT y) {
  if (n <= 0) return 0;
  OpA op;

  const int i0 = 0 * _as0, i1 = 1 * _as0, i2 = 2 * _as0;

  // unroll by rows
  ValueYType y_0 = 0, y_1 = 0, y_2 = 0;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = 0; j < n; ++j) {
    const int jj         = j * _as1;
    const ValueXType x_j = x[j * _xs0];

    y_0 += op(A[i0 + jj]) * x_j;
    y_1 += op(A[i1 + jj]) * x_j;
    y_2 += op(A[i2 + jj]) * x_j;
  }

  y[0 * _ys0] += alpha * y_0;
  y[1 * _ys0] += alpha * y_1;
  y[2 * _ys0] += alpha * y_2;

  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<2>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int n,
                                                                     ValueYType *KOKKOS_RESTRICT y) {
  if (n <= 0) return 0;
  OpA op;

  const int i0 = 0 * _as0, i1 = 1 * _as0;

  // unroll by rows
  ValueYType y_0 = 0, y_1 = 0;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = 0; j < n; ++j) {
    const int jj         = j * _as1;
    const ValueXType x_j = x[j * _xs0];

    y_0 += op(A[i0 + jj]) * x_j;
    y_1 += op(A[i1 + jj]) * x_j;
  }

  y[0 * _ys0] += alpha * y_0;
  y[1 * _ys0] += alpha * y_1;

  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<1>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int n,
                                                                     ValueYType *KOKKOS_RESTRICT y) {
  if (n <= 0) return 0;
  OpA op;

  // unroll by rows
  ValueYType y_0 = 0;

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int j = 0; j < n; ++j) y_0 += op(A[j * _as1]) * x[j * _xs0];

  y[0] += alpha * y_0;

  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<5>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int m,
                                                                     const int n, ValueYType *KOKKOS_RESTRICT y) {
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 5: {
      InnerMultipleDotProduct<5> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 4: {
      InnerMultipleDotProduct<4> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 3: {
      InnerMultipleDotProduct<3> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 2: {
      InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 1: {
      InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
  }
  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>
KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<4>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int m,
                                                                     const int n, ValueYType *KOKKOS_RESTRICT y) {
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 4: {
      InnerMultipleDotProduct<4> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 3: {
      InnerMultipleDotProduct<3> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 2: {
      InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 1: {
      InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
  }
  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>

KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<3>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int m,
                                                                     const int n, ValueYType *KOKKOS_RESTRICT y) {
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 3: {
      InnerMultipleDotProduct<3> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 2: {
      InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 1: {
      InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
  }
  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>

KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<2>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int m,
                                                                     const int n, ValueYType *KOKKOS_RESTRICT y) {
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 2: {
      InnerMultipleDotProduct<2> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
    case 1: {
      InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
  }
  return 0;
}

template <>
template <typename OpA, typename ScalarType, typename ValueAType, typename ValueXType, typename ValueYType>

KOKKOS_INLINE_FUNCTION int InnerMultipleDotProduct<1>::serial_invoke(const ScalarType alpha,
                                                                     const ValueAType *KOKKOS_RESTRICT A,
                                                                     const ValueXType *KOKKOS_RESTRICT x, const int m,
                                                                     const int n, ValueYType *KOKKOS_RESTRICT y) {
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 1: {
      InnerMultipleDotProduct<1> inner(_as0, _as1, _xs0, _ys0);
      inner.serial_invoke<OpA>(alpha, A, x, n, y);
      break;
    }
  }
  return 0;
}
}  // namespace Impl
}  // namespace KokkosBlas

#endif  // __KOKKOSBLAS_INNER_MULTIPLE_DOT_PRODUCT_SERIAL_IMPL_HPP__
