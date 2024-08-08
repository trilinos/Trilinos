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
#ifndef __KOKKOSBATCHED_INNER_TRSM_SERIAL_IMPL_HPP__
#define __KOKKOSBATCHED_INNER_TRSM_SERIAL_IMPL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_InnerTrsm_Decl.hpp"

namespace KokkosBatched {

///
/// Fixed size TRSM
/// ================
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1], a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1],
                  a_30 = A[3 * _as0 + 0 * _as1], a_31 = A[3 * _as0 + 1 * _as1], a_32 = A[3 * _as0 + 2 * _as1],
                  a_40 = A[4 * _as0 + 0 * _as1], a_41 = A[4 * _as0 + 1 * _as1], a_42 = A[4 * _as0 + 2 * _as1],
                  a_43 = A[4 * _as0 + 3 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p, ValueType &b_4p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];
    b_4p = B[4 * _bs0 + p * _bs1];

    // 0 iteration
    b_1p -= a_10 * b_0p;
    b_2p -= a_20 * b_0p;
    b_3p -= a_30 * b_0p;
    b_4p -= a_40 * b_0p;

    // 1 iteration
    b_2p -= a_21 * b_1p;
    b_3p -= a_31 * b_1p;
    b_4p -= a_41 * b_1p;

    // 2 iteration
    b_3p -= a_32 * b_2p;
    b_4p -= a_42 * b_2p;

    // 3 iteration
    b_4p -= a_43 * b_3p;

    // store
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
    B[3 * _bs0 + p * _bs1] = b_3p;
    B[4 * _bs0 + p * _bs1] = b_4p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[5];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3], b_p[4]);
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1], a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1],
                  a_30 = A[3 * _as0 + 0 * _as1], a_31 = A[3 * _as0 + 1 * _as1], a_32 = A[3 * _as0 + 2 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];

    // 0 iteration
    b_1p -= a_10 * b_0p;
    b_2p -= a_20 * b_0p;
    b_3p -= a_30 * b_0p;

    // 1 iteration
    b_2p -= a_21 * b_1p;
    b_3p -= a_31 * b_1p;

    // 2 iteration
    b_3p -= a_32 * b_2p;

    // store
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
    B[3 * _bs0 + p * _bs1] = b_3p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[4];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1], a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];

    // 0 iteration
    b_1p -= a_10 * b_0p;
    b_2p -= a_20 * b_0p;

    // 1 iteration
    b_2p -= a_21 * b_1p;

    // store
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[3];
    trsv(p, b_p[0], b_p[1], b_p[2]);
  }
  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];

    // 0 iteration
    b_1p -= a_10 * b_0p;

    // store
    B[1 * _bs0 + p * _bs1] = b_1p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[2];
    trsv(p, b_p[0], b_p[1]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT /* A */,
                                                                        const int /* n */,
                                                                        /**/ ValueType *KOKKOS_RESTRICT /* B */) {
  return 0;
}

///
/// TRSM
/// ====
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 5) Kokkos::abort("InnerTrsmLeftLowerUnitDiag<5>::serial_invoke, assert failure (m<=5)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 5: {
      InnerTrsmLeftLowerUnitDiag<5> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 4: {
      InnerTrsmLeftLowerUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftLowerUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftLowerUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 4) Kokkos::abort("InnerTrsmLeftLowerUnitDiag<4>::serial_invoke, assert failure (m<=4)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 4: {
      InnerTrsmLeftLowerUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftLowerUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftLowerUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 3) Kokkos::abort("InnerTrsmLeftLowerUnitDiag<3>::serial_invoke, assert failure (m<=3)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 3: {
      InnerTrsmLeftLowerUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftLowerUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 2) Kokkos::abort("InnerTrsmLeftLowerUnitDiag<2>::serial_invoke, assert failure (m<=2)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 2: {
      InnerTrsmLeftLowerUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 1) Kokkos::abort("InnerTrsmLeftLowerUnitDiag<1>::serial_invoke, assert failure (m<=1)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 1: {
      InnerTrsmLeftLowerUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}

///
/// Fixed size TRSM
/// ================
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1], a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1],
                  a_30 = A[3 * _as0 + 0 * _as1], a_31 = A[3 * _as0 + 1 * _as1], a_32 = A[3 * _as0 + 2 * _as1],
                  a_40 = A[4 * _as0 + 0 * _as1], a_41 = A[4 * _as0 + 1 * _as1], a_42 = A[4 * _as0 + 2 * _as1],
                  a_43 = A[4 * _as0 + 3 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1],
  //   a_22 = A[2*_as0+2*_as1],
  //   a_33 = A[3*_as0+3*_as1],
  //   a_44 = A[4*_as0+4*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1],
                  inv_a_22 = static_cast<ValueType>(1.0) / A[2 * _as0 + 2 * _as1],
                  inv_a_33 = static_cast<ValueType>(1.0) / A[3 * _as0 + 3 * _as1],
                  inv_a_44 = static_cast<ValueType>(1.0) / A[4 * _as0 + 4 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p, ValueType &b_4p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];
    b_4p = B[4 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00;*/
    b_1p -= a_10 * b_0p;
    b_2p -= a_20 * b_0p;
    b_3p -= a_30 * b_0p;
    b_4p -= a_40 * b_0p;

    // 1 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */
    b_2p -= a_21 * b_1p;
    b_3p -= a_31 * b_1p;
    b_4p -= a_41 * b_1p;

    // 2 iteration
    b_2p *= inv_a_22; /* b_2p /= a_22; */
    b_3p -= a_32 * b_2p;
    b_4p -= a_42 * b_2p;

    // 3 iteration
    b_3p *= inv_a_33; /* b_3p /= a_33; */
    b_4p -= a_43 * b_3p;

    // 4 iteration
    b_4p *= inv_a_44; /* b_4p /= a_44; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
    B[3 * _bs0 + p * _bs1] = b_3p;
    B[4 * _bs0 + p * _bs1] = b_4p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[5];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3], b_p[4]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1], a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1],
                  a_30 = A[3 * _as0 + 0 * _as1], a_31 = A[3 * _as0 + 1 * _as1], a_32 = A[3 * _as0 + 2 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1],
  //   a_22 = A[2*_as0+2*_as1],
  //   a_33 = A[3*_as0+3*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1],
                  inv_a_22 = static_cast<ValueType>(1.0) / A[2 * _as0 + 2 * _as1],
                  inv_a_33 = static_cast<ValueType>(1.0) / A[3 * _as0 + 3 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00;*/
    b_1p -= a_10 * b_0p;
    b_2p -= a_20 * b_0p;
    b_3p -= a_30 * b_0p;

    // 1 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */
    b_2p -= a_21 * b_1p;
    b_3p -= a_31 * b_1p;

    // 2 iteration
    b_2p *= inv_a_22; /* b_2p /= a_22; */
    b_3p -= a_32 * b_2p;

    // 3 iteration
    b_3p *= inv_a_33; /* b_3p /= a_33; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
    B[3 * _bs0 + p * _bs1] = b_3p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[4];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1], a_20 = A[2 * _as0 + 0 * _as1], a_21 = A[2 * _as0 + 1 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1],
  //   a_22 = A[2*_as0+2*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1],
                  inv_a_22 = static_cast<ValueType>(1.0) / A[2 * _as0 + 2 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00;*/
    b_1p -= a_10 * b_0p;
    b_2p -= a_20 * b_0p;

    // 1 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */
    b_2p -= a_21 * b_1p;

    // 2 iteration
    b_2p *= inv_a_22; /* b_2p /= a_22; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[3];
    trsv(p, b_p[0], b_p[1], b_p[2]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_10 = A[1 * _as0 + 0 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00;*/
    b_1p -= a_10 * b_0p;

    // 1 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[2];
    trsv(p, b_p[0], b_p[1]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1];

  auto trsv = [&](const int p, ValueType & /* b_0p */) {
    B[0 * _bs0 + p * _bs1] *= inv_a_00; /* b_0p /= a_00;*/
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p;
    trsv(p, b_p);
  }

  return 0;
}

///
/// TRSM
/// ==============
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 5)
    Kokkos::abort(
        "InnerTrsmLeftLowerNonUnitDiag<5>::serial_invoke, assert failure "
        "(m<=5)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 5: {
      InnerTrsmLeftLowerNonUnitDiag<5> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 4: {
      InnerTrsmLeftLowerNonUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftLowerNonUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftLowerNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 4)
    Kokkos::abort(
        "InnerTrsmLeftLowerNonUnitDiag<4>::serial_invoke, assert failure "
        "(m<=4)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 4: {
      InnerTrsmLeftLowerNonUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftLowerNonUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftLowerNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 3)
    Kokkos::abort(
        "InnerTrsmLeftLowerNonUnitDiag<3>::serial_invoke, assert failure "
        "(m<=3)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 3: {
      InnerTrsmLeftLowerNonUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftLowerNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 2)
    Kokkos::abort(
        "InnerTrsmLeftLowerNonUnitDiag<2>::serial_invoke, assert failure "
        "(m<=2)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 2: {
      InnerTrsmLeftLowerNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftLowerNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftLowerNonUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 1)
    Kokkos::abort(
        "InnerTrsmLeftLowerNonUnitDiag<1>::serial_invoke, assert failure "
        "(m<=1)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 1: {
      InnerTrsmLeftLowerNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}

///
/// Fixed size TRSM
/// ================
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1], a_03 = A[0 * _as0 + 3 * _as1],
                  a_04      = A[0 * _as0 + 4 * _as1],
                  /**/ a_12 = A[1 * _as0 + 2 * _as1], a_13 = A[1 * _as0 + 3 * _as1], a_14 = A[1 * _as0 + 4 * _as1],
                  /**/ a_23 = A[2 * _as0 + 3 * _as1], a_24 = A[2 * _as0 + 4 * _as1],
                  /**/ a_34 = A[3 * _as0 + 4 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p, ValueType &b_4p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];
    b_4p = B[4 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p -= a_04 * b_4p;
    b_1p -= a_14 * b_4p;
    b_2p -= a_24 * b_4p;
    b_3p -= a_34 * b_4p;

    // 1 iteration
    b_0p -= a_03 * b_3p;
    b_1p -= a_13 * b_3p;
    b_2p -= a_23 * b_3p;

    // 2 iteration
    b_0p -= a_02 * b_2p;
    b_1p -= a_12 * b_2p;

    // 1 iteration
    b_0p -= a_01 * b_1p;

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
    B[3 * _bs0 + p * _bs1] = b_3p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[5];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3], b_p[4]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1], a_03 = A[0 * _as0 + 3 * _as1],
                  /**/ a_12 = A[1 * _as0 + 2 * _as1], a_13 = A[1 * _as0 + 3 * _as1],
                  /**/ a_23 = A[2 * _as0 + 3 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p -= a_03 * b_3p;
    b_1p -= a_13 * b_3p;
    b_2p -= a_23 * b_3p;

    // 1 iteration
    b_0p -= a_02 * b_2p;
    b_1p -= a_12 * b_2p;

    // 2 iteration
    b_0p -= a_01 * b_1p;

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[4];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1],
                  /**/ a_12 = A[1 * _as0 + 2 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p -= a_02 * b_2p;
    b_1p -= a_12 * b_2p;

    // 1 iteration
    b_0p -= a_01 * b_1p;

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[3];
    trsv(p, b_p[0], b_p[1], b_p[2]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];

    // 0 iteration
    b_0p -= a_01 * b_1p;

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[2];
    trsv(p, b_p[0], b_p[1]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT /* A */,
                                                                        const int /* n */,
                                                                        /**/ ValueType *KOKKOS_RESTRICT /* B */) {
  return 0;
}

///
/// TRSM
/// ====
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 5) Kokkos::abort("InnerTrsmLeftUpperUnitDiag<5>::serial_invoke, assert failure (m<=5)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 5: {
      InnerTrsmLeftUpperUnitDiag<5> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 4: {
      InnerTrsmLeftUpperUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftUpperUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftUpperUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 4) Kokkos::abort("InnerTrsmLeftUpperUnitDiag<4>::serial_invoke, assert failure (m<=4)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 4: {
      InnerTrsmLeftUpperUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftUpperUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftUpperUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 3) Kokkos::abort("InnerTrsmLeftUpperUnitDiag<3>::serial_invoke, assert failure (m<=3)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 3: {
      InnerTrsmLeftUpperUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftUpperUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 2) Kokkos::abort("InnerTrsmLeftUpperUnitDiag<2>::serial_invoke, assert failure (m<=2)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 2: {
      InnerTrsmLeftUpperUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m,
                                                                        const int n,
                                                                        /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 1) Kokkos::abort("InnerTrsmLeftUpperUnitDiag<1>::serial_invoke, assert failure (m<=1)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 1: {
      InnerTrsmLeftUpperUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}

///
/// Fixed size TRSM
/// ================
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1], a_03 = A[0 * _as0 + 3 * _as1],
                  a_04      = A[0 * _as0 + 4 * _as1],
                  /**/ a_12 = A[1 * _as0 + 2 * _as1], a_13 = A[1 * _as0 + 3 * _as1], a_14 = A[1 * _as0 + 4 * _as1],
                  /**/ a_23 = A[2 * _as0 + 3 * _as1], a_24 = A[2 * _as0 + 4 * _as1],
                  /**/ a_34 = A[3 * _as0 + 4 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1],
  //   a_22 = A[2*_as0+2*_as1],
  //   a_33 = A[3*_as0+3*_as1],
  //   a_44 = A[4*_as0+4*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1],
                  inv_a_22 = static_cast<ValueType>(1.0) / A[2 * _as0 + 2 * _as1],
                  inv_a_33 = static_cast<ValueType>(1.0) / A[3 * _as0 + 3 * _as1],
                  inv_a_44 = static_cast<ValueType>(1.0) / A[4 * _as0 + 4 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p, ValueType &b_4p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];
    b_4p = B[4 * _bs0 + p * _bs1];

    // 0 iteration
    b_4p *= inv_a_44; /* b_4p /= a_44;*/
    b_3p -= a_34 * b_4p;
    b_2p -= a_24 * b_4p;
    b_1p -= a_14 * b_4p;
    b_0p -= a_04 * b_4p;

    // 1 iterationls
    b_3p *= inv_a_33; /* b_3p /= a_33;*/
    b_2p -= a_23 * b_3p;
    b_1p -= a_13 * b_3p;
    b_0p -= a_03 * b_3p;

    // 2 iteration
    b_2p *= inv_a_22; /* b_2p /= a_22; */
    b_1p -= a_12 * b_2p;
    b_0p -= a_02 * b_2p;

    // 3 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */
    b_0p -= a_01 * b_1p;

    // 4 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
    B[3 * _bs0 + p * _bs1] = b_3p;
    B[4 * _bs0 + p * _bs1] = b_4p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[5];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3], b_p[4]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1], a_03 = A[0 * _as0 + 3 * _as1],
                  /**/ a_12 = A[1 * _as0 + 2 * _as1], a_13 = A[1 * _as0 + 3 * _as1],
                  /**/ a_23 = A[2 * _as0 + 3 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1],
  //   a_22 = A[2*_as0+2*_as1],
  //   a_33 = A[3*_as0+3*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1],
                  inv_a_22 = static_cast<ValueType>(1.0) / A[2 * _as0 + 2 * _as1],
                  inv_a_33 = static_cast<ValueType>(1.0) / A[3 * _as0 + 3 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p, ValueType &b_3p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];
    b_3p = B[3 * _bs0 + p * _bs1];

    // 0 iteration
    b_3p *= inv_a_33; /* b_3p /= a_33;*/
    b_2p -= a_23 * b_3p;
    b_1p -= a_13 * b_3p;
    b_0p -= a_03 * b_3p;

    // 1 iteration
    b_2p *= inv_a_22; /* b_2p /= a_22; */
    b_1p -= a_12 * b_2p;
    b_0p -= a_02 * b_2p;

    // 2 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */
    b_0p -= a_01 * b_1p;

    // 3 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
    B[3 * _bs0 + p * _bs1] = b_3p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[4];
    trsv(p, b_p[0], b_p[1], b_p[2], b_p[3]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1], a_02 = A[0 * _as0 + 2 * _as1],
                  /**/ a_12 = A[1 * _as0 + 2 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1],
  //   a_22 = A[2*_as0+2*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1],
                  inv_a_22 = static_cast<ValueType>(1.0) / A[2 * _as0 + 2 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p, ValueType &b_2p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];
    b_2p = B[2 * _bs0 + p * _bs1];

    // 0 iteration
    b_2p *= inv_a_22; /* b_2p /= a_22; */
    b_1p -= a_12 * b_2p;
    b_0p -= a_02 * b_2p;

    // 1 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */
    b_0p -= a_01 * b_1p;

    // 2 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
    B[2 * _bs0 + p * _bs1] = b_2p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[3];
    trsv(p, b_p[0], b_p[1], b_p[2]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  const ValueType a_01 = A[0 * _as0 + 1 * _as1];

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1],
  //   a_11 = A[1*_as0+1*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1],
                  inv_a_11 = static_cast<ValueType>(1.0) / A[1 * _as0 + 1 * _as1];

  auto trsv = [&](const int p, ValueType &b_0p, ValueType &b_1p) {
    // load
    b_0p = B[0 * _bs0 + p * _bs1];
    b_1p = B[1 * _bs0 + p * _bs1];

    // 2 iteration
    b_1p *= inv_a_11; /* b_1p /= a_11; */
    b_0p -= a_01 * b_1p;

    // 3 iteration
    b_0p *= inv_a_00; /* b_0p /= a_00; */

    // store
    B[0 * _bs0 + p * _bs1] = b_0p;
    B[1 * _bs0 + p * _bs1] = b_1p;
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p[2];
    trsv(p, b_p[0], b_p[1]);
  }

  return 0;
}

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (n <= 0) return 0;

  // const ValueType
  //   a_00 = A[0*_as0+0*_as1];

  const ValueType inv_a_00 = static_cast<ValueType>(1.0) / A[0 * _as0 + 0 * _as1];

  auto trsv = [&](const int p, ValueType & /* b_0p */) {
    // 0 iteration
    B[0 * _bs0 + p * _bs1] *= inv_a_00; /* b_0p /= a_00; */
  };

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
  for (int p = 0; p < n; ++p) {
    ValueType b_p;
    trsv(p, b_p);
  }

  return 0;
}

///
/// TRSM
/// ====
/// L(m x m) X(m x n) = B (m x n)

template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<5>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 5)
    Kokkos::abort(
        "InnerTrsmLeftUpperNonUnitDiag<5>::serial_invoke, assert failure "
        "(m<=5)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 5: {
      InnerTrsmLeftUpperNonUnitDiag<5> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 4: {
      InnerTrsmLeftUpperNonUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftUpperNonUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftUpperNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<4>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 4)
    Kokkos::abort(
        "InnerTrsmLeftUpperNonUnitDiag<4>::serial_invoke, assert failure "
        "(m<=4)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 4: {
      InnerTrsmLeftUpperNonUnitDiag<4> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 3: {
      InnerTrsmLeftUpperNonUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftUpperNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<3>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 3)
    Kokkos::abort(
        "InnerTrsmLeftUpperNonUnitDiag<3>::serial_invoke, assert failure "
        "(m<=3)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 3: {
      InnerTrsmLeftUpperNonUnitDiag<3> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 2: {
      InnerTrsmLeftUpperNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<2>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 2)
    Kokkos::abort(
        "InnerTrsmLeftUpperNonUnitDiag<2>::serial_invoke, assert failure "
        "(m<=2)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 2: {
      InnerTrsmLeftUpperNonUnitDiag<2> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
    case 1: {
      InnerTrsmLeftUpperNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}
template <>
template <typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerTrsmLeftUpperNonUnitDiag<1>::serial_invoke(const ValueType *KOKKOS_RESTRICT A,
                                                                           const int m, const int n,
                                                                           /**/ ValueType *KOKKOS_RESTRICT B) {
  if (m > 1)
    Kokkos::abort(
        "InnerTrsmLeftUpperNonUnitDiag<1>::serial_invoke, assert failure "
        "(m<=1)");
  if (m <= 0 || n <= 0) return 0;
  switch (m) {
    case 1: {
      InnerTrsmLeftUpperNonUnitDiag<1> inner(_as0, _as1, _bs0, _bs1);
      inner.serial_invoke(A, n, B);
      break;
    }
  }
  return 0;
}

}  // namespace KokkosBatched

#endif
