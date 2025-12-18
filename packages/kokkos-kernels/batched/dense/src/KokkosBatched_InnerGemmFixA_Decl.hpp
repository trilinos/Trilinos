// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOSBATCHED_INNER_GEMM_FIX_A_DECL_HPP
#define KOKKOSBATCHED_INNER_GEMM_FIX_A_DECL_HPP

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {

template <int mb, int nb>
struct InnerGemmFixA {
  const int _as0, _as1, _bs0, _bs1, _cs0, _cs1;

  KOKKOS_INLINE_FUNCTION
  InnerGemmFixA(const int as0, const int as1, const int bs0, const int bs1, const int cs0, const int cs1)
      : _as0(as0), _as1(as1), _bs0(bs0), _bs1(bs1), _cs0(cs0), _cs1(cs1) {}

  // serial rank update
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const ValueType *KOKKOS_RESTRICT B, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT C);

  // serial rank update for remainder
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const ValueType *KOKKOS_RESTRICT B, const int m, const int n, const int k,
                                           /**/ ValueType *KOKKOS_RESTRICT C);
};
}  // namespace KokkosBatched

#endif
