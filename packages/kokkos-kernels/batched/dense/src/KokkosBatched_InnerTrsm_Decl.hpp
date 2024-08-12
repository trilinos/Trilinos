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
#ifndef __KOKKOSBATCHED_INNER_TRSM_DECL_HPP__
#define __KOKKOSBATCHED_INNER_TRSM_DECL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace KokkosBatched {

// specialized for different m and n
// Solve L(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftLowerUnitDiag {
  const int _as0, _as1, _bs0, _bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftLowerUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : _as0(as0), _as1(as1), _bs0(bs0), _bs1(bs1) {}

  // trisolve
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

// specialized for different m and n
// Solve L(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftLowerNonUnitDiag {
  const int _as0, _as1, _bs0, _bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftLowerNonUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : _as0(as0), _as1(as1), _bs0(bs0), _bs1(bs1) {}

  // trisolve
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

// specialized for different m and n
// Solve U(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftUpperUnitDiag {
  const int _as0, _as1, _bs0, _bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftUpperUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : _as0(as0), _as1(as1), _bs0(bs0), _bs1(bs1) {}

  // trisolve
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

// specialized for different m and n
// Solve U(m x m) X(m x n) = B(m x n)
template <int bmn>
struct InnerTrsmLeftUpperNonUnitDiag {
  const int _as0, _as1, _bs0, _bs1;

  KOKKOS_INLINE_FUNCTION
  InnerTrsmLeftUpperNonUnitDiag(const int as0, const int as1, const int bs0, const int bs1)
      : _as0(as0), _as1(as1), _bs0(bs0), _bs1(bs1) {}

  // trisolve
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);

  // for remainder
  template <typename ValueType>
  KOKKOS_INLINE_FUNCTION int serial_invoke(const ValueType *KOKKOS_RESTRICT A, const int m, const int n,
                                           /**/ ValueType *KOKKOS_RESTRICT B);
};

}  // namespace KokkosBatched

#endif
