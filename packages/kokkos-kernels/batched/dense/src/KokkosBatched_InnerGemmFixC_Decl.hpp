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
#ifndef __KOKKOSBATCHED_INNER_GEMM_FIX_C_DECL_HPP__
#define __KOKKOSBATCHED_INNER_GEMM_FIX_C_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {

template <int mb = 0, int nb = 0>
struct InnerGemmFixC {
  const int _as0, _as1, _bs0, _bs1, _cs0, _cs1;

  KOKKOS_INLINE_FUNCTION
  InnerGemmFixC(const int as0, const int as1, const int bs0, const int bs1, const int cs0, const int cs1)
      : _as0(as0), _as1(as1), _bs0(bs0), _bs1(bs1), _cs0(cs0), _cs1(cs1) {}

  // serial rank update
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const ValueType *KOKKOS_RESTRICT B, const int k,
                                           /**/ ValueType *KOKKOS_RESTRICT C);

  // serial rank update for remainder
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const ValueType *KOKKOS_RESTRICT B, const int m, const int k,
                                           /**/ ValueType *KOKKOS_RESTRICT C);

  // serial rank update for remainder
  template <typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ScalarType alpha, const ValueType *KOKKOS_RESTRICT A,
                                           const ValueType *KOKKOS_RESTRICT B, const int m, const int n, const int k,
                                           /**/ ValueType *KOKKOS_RESTRICT C);

  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int team_invoke(const MemberType &member, const ScalarType alpha,
                                         const ValueType *KOKKOS_RESTRICT A, const ValueType *KOKKOS_RESTRICT B,
                                         const int k,
                                         /**/ ValueType *KOKKOS_RESTRICT C);

  // team rank update for remainder
  template <typename MemberType, typename ScalarType, typename ValueType>
  KOKKOS_INLINE_FUNCTION int team_invoke(const MemberType &member, const ScalarType alpha,
                                         const ValueType *KOKKOS_RESTRICT A, const ValueType *KOKKOS_RESTRICT B,
                                         const int m, const int n, const int k,
                                         /**/ ValueType *KOKKOS_RESTRICT C);
};
}  // namespace KokkosBatched

#endif
