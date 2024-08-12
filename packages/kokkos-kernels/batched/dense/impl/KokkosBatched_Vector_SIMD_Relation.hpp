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
#ifndef __KOKKOSBATCHED_VECTOR_SIMD_RELATION_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_RELATION_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Complex.hpp"

namespace KokkosBatched {

// vector, vector

#undef KOKKOSBATCHED_RELATION_OPERATOR
#define KOKKOSBATCHED_RELATION_OPERATOR(op)                                                      \
  template <typename T1, typename T2, int l>                                                     \
  KOKKOS_INLINE_FUNCTION const Vector<SIMD<bool>, l> operator op(const Vector<SIMD<T1>, l> &a,   \
                                                                 const Vector<SIMD<T2>, l> &b) { \
    Vector<SIMD<bool>, l> r_val;                                                                 \
    for (int i = 0; i < l; ++i) r_val[i] = a[i] op b[i];                                         \
    return r_val;                                                                                \
  }

KOKKOSBATCHED_RELATION_OPERATOR(<)
KOKKOSBATCHED_RELATION_OPERATOR(>)
KOKKOSBATCHED_RELATION_OPERATOR(<=)
KOKKOSBATCHED_RELATION_OPERATOR(>=)
KOKKOSBATCHED_RELATION_OPERATOR(==)
KOKKOSBATCHED_RELATION_OPERATOR(!=)

// vector, scalar
#undef KOKKOSBATCHED_RELATION_OPERATOR
#define KOKKOSBATCHED_RELATION_OPERATOR(op)                                                                   \
  template <typename T1, typename T2, int l>                                                                  \
  KOKKOS_INLINE_FUNCTION const Vector<SIMD<bool>, l> operator op(const Vector<SIMD<T1>, l> &a, const T2 &b) { \
    Vector<SIMD<bool>, l> r_val;                                                                              \
    for (int i = 0; i < l; ++i) r_val[i] = a[i] op b;                                                         \
    return r_val;                                                                                             \
  }

KOKKOSBATCHED_RELATION_OPERATOR(<)
KOKKOSBATCHED_RELATION_OPERATOR(>)
KOKKOSBATCHED_RELATION_OPERATOR(<=)
KOKKOSBATCHED_RELATION_OPERATOR(>=)
KOKKOSBATCHED_RELATION_OPERATOR(==)
KOKKOSBATCHED_RELATION_OPERATOR(!=)

// scalar, vector
#undef KOKKOSBATCHED_RELATION_OPERATOR
#define KOKKOSBATCHED_RELATION_OPERATOR(op)                                                                   \
  template <typename T1, typename T2, int l>                                                                  \
  KOKKOS_INLINE_FUNCTION const Vector<SIMD<bool>, l> operator op(const T1 &a, const Vector<SIMD<T2>, l> &b) { \
    Vector<SIMD<bool>, l> r_val;                                                                              \
    for (int i = 0; i < l; ++i) r_val[i] = a op b[i];                                                         \
    return r_val;                                                                                             \
  }

KOKKOSBATCHED_RELATION_OPERATOR(<)
KOKKOSBATCHED_RELATION_OPERATOR(>)
KOKKOSBATCHED_RELATION_OPERATOR(<=)
KOKKOSBATCHED_RELATION_OPERATOR(>=)
KOKKOSBATCHED_RELATION_OPERATOR(==)
KOKKOSBATCHED_RELATION_OPERATOR(!=)

#undef KOKKOSBATCHED_RELATION_OPERATOR
}  // namespace KokkosBatched

#endif
