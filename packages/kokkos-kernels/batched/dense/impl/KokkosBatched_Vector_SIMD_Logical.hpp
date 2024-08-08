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
#ifndef __KOKKOSBATCHED_VECTOR_SIMD_LOGICAL_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_LOGICAL_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Complex.hpp"

namespace KokkosBatched {

#define KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE(T0, T1, l)                        \
  typename std::enable_if<std::is_integral<T0>::value && std::is_integral<T1>::value, \
                          const Vector<SIMD<bool>, l> >::type

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static typename std::enable_if<std::is_integral<T>::value, const Vector<SIMD<bool>, l> >::type
operator!(const Vector<SIMD<T>, l> &a) {
  Vector<SIMD<bool>, l> r_val;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_VECTOR)
#pragma vector always
#endif
  for (int i = 0; i < l; ++i) r_val[i] = !a[i];
  return r_val;
}

template <typename T0, typename T1, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE(T0, T1, l) operator||(
    const Vector<SIMD<T0>, l> &a, const Vector<SIMD<T1>, l> &b) {
  Vector<SIMD<bool>, l> r_val;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_VECTOR)
#pragma vector always
#endif
  for (int i = 0; i < l; ++i) r_val[i] = a[i] || b[i];
  return r_val;
}

template <typename T0, typename T1, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE(T0, T1, l) operator&&(
    const Vector<SIMD<T0>, l> &a, const Vector<SIMD<T1>, l> &b) {
  Vector<SIMD<bool>, l> r_val;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_VECTOR)
#pragma vector always
#endif
  for (int i = 0; i < l; ++i) r_val[i] = a[i] && b[i];
  return r_val;
}

template <typename T0, typename T1, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE(T0, T1, l) operator||(
    const Vector<SIMD<T0>, l> &a, const T1 &b) {
  Vector<SIMD<bool>, l> r_val;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_VECTOR)
#pragma vector always
#endif
  for (int i = 0; i < l; ++i) r_val[i] = a[i] || b;
  return r_val;
}

template <typename T0, typename T1, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE(T0, T1, l) operator&&(
    const Vector<SIMD<T0>, l> &a, const T1 &b) {
  Vector<SIMD<bool>, l> r_val;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_VECTOR)
#pragma vector always
#endif
  for (int i = 0; i < l; ++i) r_val[i] = a[i] && b;
  return r_val;
}

template <typename T0, typename T1, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE(T0, T1, l) operator||(
    const T0 &a, const Vector<SIMD<T1>, l> &b) {
  Vector<SIMD<bool>, l> r_val;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_VECTOR)
#pragma vector always
#endif
  for (int i = 0; i < l; ++i) r_val[i] = a || b[i];
  return r_val;
}

template <typename T0, typename T1, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE(T0, T1, l) operator&&(
    const T0 &a, const Vector<SIMD<T1>, l> &b) {
  Vector<SIMD<bool>, l> r_val;
#if defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
#if defined(KOKKOS_ENABLE_PRAGMA_VECTOR)
#pragma vector always
#endif
  for (int i = 0; i < l; ++i) r_val[i] = a && b[i];
  return r_val;
}
#undef KOKKOSKERNELS_SIMD_LOGICAL_RETURN_BOOL_TYPE
}  // namespace KokkosBatched

#endif
