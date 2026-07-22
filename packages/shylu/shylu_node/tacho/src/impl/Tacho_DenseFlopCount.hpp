// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_DENSE_FLOP_HPP__
#define __TACHO_DENSE_FLOP_HPP__

#include "Tacho_Util.hpp"

/// \file Tacho_DenseFlopCount.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

//  FLOP counting - From LAPACK working note #41

namespace Tacho {

#define FLOP_MUL(isComplex) ((isComplex) ? (6.0) : (1.0))
#define FLOP_ADD(isComplex) ((isComplex) ? (2.0) : (1.0))

template <typename ValueType> class DenseFlopCount {
public:
  static KOKKOS_INLINE_FUNCTION double Gemm(int mm, int nn, int kk) {
    double m = (double)mm;
    double n = (double)nn;
    double k = (double)kk;
    return (FLOP_MUL(ArithTraits<ValueType>::is_complex) * (m * n * k) +
            FLOP_ADD(ArithTraits<ValueType>::is_complex) * (m * n * k));
  }

  static KOKKOS_INLINE_FUNCTION double Syrk(int kk, int nn) {
    double k = (double)kk;
    double n = (double)nn;
    return (FLOP_MUL(ArithTraits<ValueType>::is_complex) * (0.5 * k * n * (n + 1.0)) +
            FLOP_ADD(ArithTraits<ValueType>::is_complex) * (0.5 * k * n * (n + 1.0)));
  }

  static KOKKOS_INLINE_FUNCTION double TrsmLower(int mm, int nn) {
    double m = (double)mm;
    double n = (double)nn;
    return (FLOP_MUL(ArithTraits<ValueType>::is_complex) * (0.5 * n * m * (m + 1.0)) +
            FLOP_ADD(ArithTraits<ValueType>::is_complex) * (0.5 * n * m * (m - 1.0)));
  }

  static KOKKOS_INLINE_FUNCTION double TrsmUpper(int mm, int nn) {
    double m = (double)mm;
    double n = (double)nn;
    return (FLOP_MUL(ArithTraits<ValueType>::is_complex) * (0.5 * m * n * (n + 1.0)) +
            FLOP_ADD(ArithTraits<ValueType>::is_complex) * (0.5 * m * n * (n - 1.0)));
  }

  static KOKKOS_INLINE_FUNCTION double Trsm(int is_lower, int mm, int nn) {
    return (is_lower ? TrsmLower(mm, nn) : TrsmUpper(mm, nn));
  }

  static KOKKOS_INLINE_FUNCTION double LU(int mm, int nn) {
    double m = (double)mm;
    double n = (double)nn;
    if (m > n)
      return (FLOP_MUL(ArithTraits<ValueType>::is_complex) *
                  (0.5 * m * n * n - (1.0 / 6.0) * n * n * n + 0.5 * m * n - 0.5 * n * n + (2.0 / 3.0) * n) +
              FLOP_ADD(ArithTraits<ValueType>::is_complex) *
                  (0.5 * m * n * n - (1.0 / 6.0) * n * n * n - 0.5 * m * n + (1.0 / 6.0) * n));
    else
      return (FLOP_MUL(ArithTraits<ValueType>::is_complex) *
                  (0.5 * n * m * m - (1.0 / 6.0) * m * m * m + 0.5 * n * m - 0.5 * m * m + (2.0 / 3.0) * m) +
              FLOP_ADD(ArithTraits<ValueType>::is_complex) *
                  (0.5 * n * m * m - (1.0 / 6.0) * m * m * m - 0.5 * n * m + (1.0 / 6.0) * m));
  }

  static KOKKOS_INLINE_FUNCTION double Chol(int nn) {
    double n = (double)nn;
    return (FLOP_MUL(ArithTraits<ValueType>::is_complex) * ((1.0 / 6.0) * n * n * n + 0.5 * n * n + (1.0 / 3.0) * n) +
            FLOP_ADD(ArithTraits<ValueType>::is_complex) * ((1.0 / 6.0) * n * n * n - (1.0 / 6.0) * n));
  }

  static KOKKOS_INLINE_FUNCTION double LDL(int nn) {
    double n = (double)nn;
    return (FLOP_MUL(ArithTraits<ValueType>::is_complex) * ((1.0 / 3.0) * n * n * n + (2.0 / 3.0) * n) +
            FLOP_ADD(ArithTraits<ValueType>::is_complex) * ((1.0 / 3.0) * n * n * n - (1.0 / 3.0) * n));
  }
};

} // namespace Tacho

#endif
