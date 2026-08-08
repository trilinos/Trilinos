// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_INNER_GEMM_FIX_A_SERIAL_IMPL_HPP
#define KOKKOSBATCHED_INNER_GEMM_FIX_A_SERIAL_IMPL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_InnerGemmFixA_Decl.hpp"

namespace KokkosBatched {

///
/// Inner kernel (5x5)
/// ==================

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<5, 5>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_04 = opA(A[0 * m_as0 + 4 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]),
                  a_13 = opA(A[1 * m_as0 + 3 * m_as1]), a_14 = opA(A[1 * m_as0 + 4 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_23 = opA(A[2 * m_as0 + 3 * m_as1]),
                  a_24 = opA(A[2 * m_as0 + 4 * m_as1]), a_30 = opA(A[3 * m_as0 + 0 * m_as1]),
                  a_31 = opA(A[3 * m_as0 + 1 * m_as1]), a_32 = opA(A[3 * m_as0 + 2 * m_as1]),
                  a_33 = opA(A[3 * m_as0 + 3 * m_as1]), a_34 = opA(A[3 * m_as0 + 4 * m_as1]),
                  a_40 = opA(A[4 * m_as0 + 0 * m_as1]), a_41 = opA(A[4 * m_as0 + 1 * m_as1]),
                  a_42 = opA(A[4 * m_as0 + 2 * m_as1]), a_43 = opA(A[4 * m_as0 + 3 * m_as1]),
                  a_44 = opA(A[4 * m_as0 + 4 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p, b_3p, c_3p, b_4p, c_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ib4 = 4 * m_bs0, ic0 = 0 * m_cs0,
            ic1 = 1 * m_cs0, ic2 = 2 * m_cs0, ic3 = 3 * m_cs0, ic4 = 4 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);
    b_4p = opB(B[ib4 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_0p += a_04 * b_4p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;
    c_1p += a_14 * b_4p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_2p += a_23 * b_3p;
    c_2p += a_24 * b_4p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;
    c_3p += a_32 * b_2p;
    c_3p += a_33 * b_3p;
    c_3p += a_34 * b_4p;
    c_4p = a_40 * b_0p;
    c_4p += a_41 * b_1p;
    c_4p += a_42 * b_2p;
    c_4p += a_43 * b_3p;
    c_4p += a_44 * b_4p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
    C[ic4 + p * m_cs1] += alpha * c_4p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<5, 4>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]),
                  a_12 = opA(A[1 * m_as0 + 2 * m_as1]), a_13 = opA(A[1 * m_as0 + 3 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_23 = opA(A[2 * m_as0 + 3 * m_as1]),
                  a_30 = opA(A[3 * m_as0 + 0 * m_as1]), a_31 = opA(A[3 * m_as0 + 1 * m_as1]),
                  a_32 = opA(A[3 * m_as0 + 2 * m_as1]), a_33 = opA(A[3 * m_as0 + 3 * m_as1]),
                  a_40 = opA(A[4 * m_as0 + 0 * m_as1]), a_41 = opA(A[4 * m_as0 + 1 * m_as1]),
                  a_42 = opA(A[4 * m_as0 + 2 * m_as1]), a_43 = opA(A[4 * m_as0 + 3 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p, b_3p, c_3p,
      /**/ c_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0,
            ic2 = 2 * m_cs0, ic3 = 3 * m_cs0, ic4 = 4 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_2p += a_23 * b_3p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;
    c_3p += a_32 * b_2p;
    c_3p += a_33 * b_3p;
    c_4p = a_40 * b_0p;
    c_4p += a_41 * b_1p;
    c_4p += a_42 * b_2p;
    c_4p += a_43 * b_3p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
    C[ic4 + p * m_cs1] += alpha * c_4p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<5, 3>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_30 = opA(A[3 * m_as0 + 0 * m_as1]),
                  a_31 = opA(A[3 * m_as0 + 1 * m_as1]), a_32 = opA(A[3 * m_as0 + 2 * m_as1]),
                  a_40 = opA(A[4 * m_as0 + 0 * m_as1]), a_41 = opA(A[4 * m_as0 + 1 * m_as1]),
                  a_42 = opA(A[4 * m_as0 + 2 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p,
      /**/ c_3p,
      /**/ c_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0,
            ic3 = 3 * m_cs0, ic4 = 4 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;
    c_3p += a_32 * b_2p;
    c_4p = a_40 * b_0p;
    c_4p += a_41 * b_1p;
    c_4p += a_42 * b_2p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
    C[ic4 + p * m_cs1] += alpha * c_4p;
  }

  return 0;
}
template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<5, 2>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_30 = opA(A[3 * m_as0 + 0 * m_as1]), a_31 = opA(A[3 * m_as0 + 1 * m_as1]),
                  a_40 = opA(A[4 * m_as0 + 0 * m_as1]), a_41 = opA(A[4 * m_as0 + 1 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p,
      /**/ c_2p,
      /**/ c_3p,
      /**/ c_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0, ic3 = 3 * m_cs0,
            ic4 = 4 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;
    c_4p = a_40 * b_0p;
    c_4p += a_41 * b_1p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
    C[ic4 + p * m_cs1] += alpha * c_4p;
  }

  return 0;
}
template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<5, 1>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_30 = opA(A[3 * m_as0 + 0 * m_as1]),
                  a_40 = opA(A[4 * m_as0 + 0 * m_as1]);

  ValueType b_0p, c_0p,
      /**/ c_1p,
      /**/ c_2p,
      /**/ c_3p,
      /**/ c_4p;

  const int ib0 = 0 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0, ic3 = 3 * m_cs0, ic4 = 4 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_1p = a_10 * b_0p;
    c_2p = a_20 * b_0p;
    c_3p = a_30 * b_0p;
    c_4p = a_40 * b_0p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
    C[ic4 + p * m_cs1] += alpha * c_4p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<4, 5>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_04 = opA(A[0 * m_as0 + 4 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]),
                  a_13 = opA(A[1 * m_as0 + 3 * m_as1]), a_14 = opA(A[1 * m_as0 + 4 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_23 = opA(A[2 * m_as0 + 3 * m_as1]),
                  a_24 = opA(A[2 * m_as0 + 4 * m_as1]), a_30 = opA(A[3 * m_as0 + 0 * m_as1]),
                  a_31 = opA(A[3 * m_as0 + 1 * m_as1]), a_32 = opA(A[3 * m_as0 + 2 * m_as1]),
                  a_33 = opA(A[3 * m_as0 + 3 * m_as1]), a_34 = opA(A[3 * m_as0 + 4 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p, b_3p, c_3p, b_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ib4 = 4 * m_bs0, ic0 = 0 * m_cs0,
            ic1 = 1 * m_cs0, ic2 = 2 * m_cs0, ic3 = 3 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);
    b_4p = opB(B[ib4 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_0p += a_04 * b_4p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;
    c_1p += a_14 * b_4p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_2p += a_23 * b_3p;
    c_2p += a_24 * b_4p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;
    c_3p += a_32 * b_2p;
    c_3p += a_33 * b_3p;
    c_3p += a_34 * b_4p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<3, 5>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_04 = opA(A[0 * m_as0 + 4 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]),
                  a_13 = opA(A[1 * m_as0 + 3 * m_as1]), a_14 = opA(A[1 * m_as0 + 4 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_23 = opA(A[2 * m_as0 + 3 * m_as1]),
                  a_24 = opA(A[2 * m_as0 + 4 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p, b_3p, b_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ib4 = 4 * m_bs0, ic0 = 0 * m_cs0,
            ic1 = 1 * m_cs0, ic2 = 2 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);
    b_4p = opB(B[ib4 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_0p += a_04 * b_4p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;
    c_1p += a_14 * b_4p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_2p += a_23 * b_3p;
    c_2p += a_24 * b_4p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
  }

  return 0;
}
template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<2, 5>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_04 = opA(A[0 * m_as0 + 4 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]),
                  a_13 = opA(A[1 * m_as0 + 3 * m_as1]), a_14 = opA(A[1 * m_as0 + 4 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, b_3p, b_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ib4 = 4 * m_bs0, ic0 = 0 * m_cs0,
            ic1 = 1 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);
    b_4p = opB(B[ib4 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_0p += a_04 * b_4p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;
    c_1p += a_14 * b_4p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
  }

  return 0;
}
template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<1, 5>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_04 = opA(A[0 * m_as0 + 4 * m_as1]);

  ValueType b_0p, c_0p, b_1p, b_2p, b_3p, b_4p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ib4 = 4 * m_bs0, ic0 = 0 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);
    b_4p = opB(B[ib4 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_0p += a_04 * b_4p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<5, 5>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int m,
                                                              const int n, const int k,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (m <= 0 || n <= 0 || k <= 0) return 0;

  switch (m * 10 + k) {
    case 54: {
      InnerGemmFixA<5, 4> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 53: {
      InnerGemmFixA<5, 3> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 52: {
      InnerGemmFixA<5, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 51: {
      InnerGemmFixA<5, 1> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 45: {
      InnerGemmFixA<4, 5> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 35: {
      InnerGemmFixA<3, 5> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 25: {
      InnerGemmFixA<2, 5> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 15: {
      InnerGemmFixA<1, 5> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    default: {
      if (m < 5 && n < 5) {
        InnerGemmFixA<2, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
        for (int i = 0; i < m; i += 2)
          for (int p = 0; p < k; p += 2)
            inner.serial_invoke(opA, opB, alpha, A + i * m_as0 + p * m_as1, B + p * m_bs0, (i + 2 > m ? 1 : 2), n,
                                (p + 2 > k ? 1 : 2), C + i * m_cs0);
      } else {
        Kokkos::abort("InnerGemmFixA<5,5>::serial_invoke, assert failure (m<5 && n<5)");
      }
      break;
    }
  }

  return 0;
}

///
/// Inner kernel (4x4)
/// ==================

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<4, 4>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]),
                  a_12 = opA(A[1 * m_as0 + 2 * m_as1]), a_13 = opA(A[1 * m_as0 + 3 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_23 = opA(A[2 * m_as0 + 3 * m_as1]),
                  a_30 = opA(A[3 * m_as0 + 0 * m_as1]), a_31 = opA(A[3 * m_as0 + 1 * m_as1]),
                  a_32 = opA(A[3 * m_as0 + 2 * m_as1]), a_33 = opA(A[3 * m_as0 + 3 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p, b_3p, c_3p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0,
            ic2 = 2 * m_cs0, ic3 = 3 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_2p += a_23 * b_3p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;
    c_3p += a_32 * b_2p;
    c_3p += a_33 * b_3p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<4, 3>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_30 = opA(A[3 * m_as0 + 0 * m_as1]),
                  a_31 = opA(A[3 * m_as0 + 1 * m_as1]), a_32 = opA(A[3 * m_as0 + 2 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p,
      /**/ c_3p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0,
            ic3 = 3 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;
    c_3p += a_32 * b_2p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<4, 2>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_30 = opA(A[3 * m_as0 + 0 * m_as1]), a_31 = opA(A[3 * m_as0 + 1 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p,
      /**/ c_2p,
      /**/ c_3p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0, ic3 = 3 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_3p = a_30 * b_0p;
    c_3p += a_31 * b_1p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<4, 1>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_30 = opA(A[3 * m_as0 + 0 * m_as1]);

  ValueType b_0p, c_0p,
      /**/ c_1p,
      /**/ c_2p,
      /**/ c_3p;

  const int ib0 = 0 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0, ic3 = 3 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_1p = a_10 * b_0p;
    c_2p = a_20 * b_0p;
    c_3p = a_30 * b_0p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
    C[ic3 + p * m_cs1] += alpha * c_3p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<3, 4>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]),
                  a_12 = opA(A[1 * m_as0 + 2 * m_as1]), a_13 = opA(A[1 * m_as0 + 3 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]), a_23 = opA(A[2 * m_as0 + 3 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p, b_3p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0,
            ic2 = 2 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;
    c_2p += a_23 * b_3p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<2, 4>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]),
                  a_12 = opA(A[1 * m_as0 + 2 * m_as1]), a_13 = opA(A[1 * m_as0 + 3 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, b_3p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_1p += a_13 * b_3p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<1, 4>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_03 = opA(A[0 * m_as0 + 3 * m_as1]);

  ValueType b_0p, c_0p, b_1p, b_2p, b_3p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ib3 = 3 * m_bs0, ic0 = 0 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);
    b_3p = opB(B[ib3 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_0p += a_03 * b_3p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<4, 4>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int m,
                                                              const int n, const int k,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (m <= 0 || n <= 0 || k <= 0) return 0;

  switch (m * 10 + k) {
    case 44: {
      InnerGemmFixA<4, 4> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 43: {
      InnerGemmFixA<4, 3> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 42: {
      InnerGemmFixA<4, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 41: {
      InnerGemmFixA<4, 1> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 34: {
      InnerGemmFixA<3, 4> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 24: {
      InnerGemmFixA<2, 4> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 14: {
      InnerGemmFixA<1, 4> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    default: {
      if (m < 4 && n < 4) {
        InnerGemmFixA<2, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
        for (int i = 0; i < m; i += 2)
          for (int p = 0; p < k; p += 2)
            inner.serial_invoke(opA, opB, alpha, A + i * m_as0 + p * m_as1, B + p * m_bs0, (i + 2 > m ? 1 : 2), n,
                                (p + 2 > k ? 1 : 2), C + i * m_cs0);
      } else {
        Kokkos::abort("InnerGemmFixA<4,4>::serial_invoke, assert failure (m<4 && n<4)");
      }
      break;
    }
  }

  return 0;
}

///
/// Inner kernel (3x3)
/// ==================

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<3, 3>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]),
                  a_22 = opA(A[2 * m_as0 + 2 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p, c_2p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;
    c_2p += a_22 * b_2p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<3, 2>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]), a_21 = opA(A[2 * m_as0 + 1 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p,
      /**/ c_2p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_2p = a_20 * b_0p;
    c_2p += a_21 * b_1p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<3, 1>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_20 = opA(A[2 * m_as0 + 0 * m_as1]);

  ValueType b_0p, c_0p,
      /**/ c_1p,
      /**/ c_2p;

  const int ib0 = 0 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0, ic2 = 2 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_1p = a_10 * b_0p;
    c_2p = a_20 * b_0p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
    C[ic2 + p * m_cs1] += alpha * c_2p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<2, 3>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]),
                  a_11 = opA(A[1 * m_as0 + 1 * m_as1]), a_12 = opA(A[1 * m_as0 + 2 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p, b_2p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;
    c_1p += a_12 * b_2p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
  }

  return 0;
}
template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<1, 3>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_02 = opA(A[0 * m_as0 + 2 * m_as1]);

  ValueType b_0p, c_0p, b_1p, b_2p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ib2 = 2 * m_bs0, ic0 = 0 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);
    b_2p = opB(B[ib2 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_0p += a_02 * b_2p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<3, 3>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int m,
                                                              const int n, const int k,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (m <= 0 || n <= 0 || k <= 0) return 0;

  switch (m * 10 + k) {
    case 33: {
      InnerGemmFixA<3, 3> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 32: {
      InnerGemmFixA<3, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 31: {
      InnerGemmFixA<3, 1> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 23: {
      InnerGemmFixA<2, 3> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 13: {
      InnerGemmFixA<1, 3> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    default: {
      if (m < 3 && n < 3) {
        InnerGemmFixA<2, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
        for (int i = 0; i < m; i += 2)
          for (int p = 0; p < k; p += 2)
            inner.serial_invoke(opA, opB, alpha, A + i * m_as0 + p * m_as1, B + p * m_bs0, (i + 2 > m ? 1 : 2), n,
                                (p + 2 > k ? 1 : 2), C + i * m_cs0);
      } else {
        Kokkos::abort("InnerGemmFixA<3,3>::serial_invoke, assert failure (m<3 && n<3)");
      }
      break;
    }
  }

  return 0;
}

///
/// Inner kernel (2x2)
/// ==================

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<2, 2>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]),
                  a_10 = opA(A[1 * m_as0 + 0 * m_as1]), a_11 = opA(A[1 * m_as0 + 1 * m_as1]);

  ValueType b_0p, c_0p, b_1p, c_1p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;
    c_1p = a_10 * b_0p;
    c_1p += a_11 * b_1p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<2, 1>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_10 = opA(A[1 * m_as0 + 0 * m_as1]);

  ValueType b_0p, c_0p,
      /**/ c_1p;

  const int ib0 = 0 * m_bs0, ic0 = 0 * m_cs0, ic1 = 1 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_1p = a_10 * b_0p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
    C[ic1 + p * m_cs1] += alpha * c_1p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<1, 2>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]), a_01 = opA(A[0 * m_as0 + 1 * m_as1]);

  ValueType b_0p, c_0p, b_1p;

  const int ib0 = 0 * m_bs0, ib1 = 1 * m_bs0, ic0 = 0 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);
    b_1p = opB(B[ib1 + p * m_bs1]);

    c_0p = a_00 * b_0p;
    c_0p += a_01 * b_1p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
  }

  return 0;
}

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<2, 2>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int m,
                                                              const int n, const int k,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (m <= 0 || n <= 0 || k <= 0) return 0;

  switch (m * 10 + k) {
    case 22: {
      InnerGemmFixA<2, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 21: {
      InnerGemmFixA<2, 1> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 12: {
      InnerGemmFixA<1, 2> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    case 11: {
      InnerGemmFixA<1, 1> inner(m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1);
      inner.serial_invoke(opA, opB, alpha, A, B, n, C);
      break;
    }
    default: {
      Kokkos::abort("InnerGemmFixA<2,2>::serial_invoke, assert failure (m<2 && n<2)");
      break;
    }
  }

  return 0;
}

///
/// Inner kernel (1x1)
/// ==================

template <>
template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
KOKKOS_INLINE_FUNCTION int InnerGemmFixA<1, 1>::serial_invoke(OpA opA, OpB opB, const ScalarType alpha,
                                                              const ValueType *KOKKOS_RESTRICT A,
                                                              const ValueType *KOKKOS_RESTRICT B, const int n,
                                                              /**/ ValueType *KOKKOS_RESTRICT C) {
  if (n <= 0) return 0;

  const ValueType a_00 = opA(A[0 * m_as0 + 0 * m_as1]);

  ValueType b_0p, c_0p;

  const int ib0 = 0 * m_bs0, ic0 = 0 * m_cs0;

  for (int p = 0; p < n; ++p) {
    b_0p = opB(B[ib0 + p * m_bs1]);

    c_0p = a_00 * b_0p;

    C[ic0 + p * m_cs1] += alpha * c_0p;
  }

  return 0;
}

}  // namespace KokkosBatched

#endif
