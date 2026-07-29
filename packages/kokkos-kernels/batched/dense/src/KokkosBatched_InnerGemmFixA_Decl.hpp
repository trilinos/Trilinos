// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_INNER_GEMM_FIX_A_DECL_HPP
#define KOKKOSBATCHED_INNER_GEMM_FIX_A_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {

template <int mb, int nb>
struct InnerGemmFixA {
  const int m_as0, m_as1, m_bs0, m_bs1, m_cs0, m_cs1;

  KOKKOS_INLINE_FUNCTION
  InnerGemmFixA(const int as0, const int as1, const int bs0, const int bs1, const int cs0, const int cs1)
      : m_as0(as0), m_as1(as1), m_bs0(bs0), m_bs1(bs1), m_cs0(cs0), m_cs1(cs1) {}

  // serial rank update
  template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(OpA opA, OpB opB, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const ValueType *KOKKOS_RESTRICT B, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT C);

  // serial rank update for remainder
  template <typename OpA, typename OpB, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(OpA opA, OpB opB, const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const ValueType *KOKKOS_RESTRICT B, const int m, const int n, const int k,
                                           /**/ ValueType *KOKKOS_RESTRICT C);
};
}  // namespace KokkosBatched

#endif
