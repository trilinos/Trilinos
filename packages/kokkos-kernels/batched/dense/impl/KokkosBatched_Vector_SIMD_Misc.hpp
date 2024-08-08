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
#ifndef __KOKKOSBATCHED_VECTOR_SIMD_MISC_HPP__
#define __KOKKOSBATCHED_VECTOR_SIMD_MISC_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Kokkos_Complex.hpp"

namespace KokkosBatched {

#define KOKKOSKERNELS_SIMD_MISC_RETURN_TYPE(T, l) Vector<SIMD<T>, l>
#define KOKKOSKERNELS_SIMD_MISC_CONVERTIBLE_RETURN_VOID_TYPE(T0, T1, T2, l) void
// typename std::enable_if<std::is_convertible< T1 , T0 >::value &&
// std::is_convertible< T2 , T0 >::value, void >::type

// scalar, scalar

template <typename T>
KOKKOS_INLINE_FUNCTION static T conditional_assign(const bool cond, const T &if_true_val, const T &if_false_val) {
  return cond ? if_true_val : if_false_val;
}

template <typename T0, typename T1, typename T2>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_MISC_CONVERTIBLE_RETURN_VOID_TYPE(T0, T1, T2, l)
    conditional_assign(/* */ T0 &r_val, const bool cond, const T1 &if_true_val, const T2 &if_false_val) {
  r_val = cond ? if_true_val : if_false_val;
}

// vector, scalar

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_MISC_RETURN_TYPE(T, l)
    conditional_assign(const Vector<SIMD<bool>, l> &cond, const Vector<SIMD<T>, l> &if_true_val,
                       const T &if_false_val) {
  Vector<SIMD<T>, l> r_val;
  for (int i = 0; i < l; ++i) r_val[i] = cond[i] ? if_true_val[i] : if_false_val;
  return r_val;
}

template <typename T0, typename T1, typename T2, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_MISC_CONVERTIBLE_RETURN_VOID_TYPE(T0, T1, T2, l)
    conditional_assign(/* */ Vector<SIMD<T0>, l> &r_val, const Vector<SIMD<bool>, l> &cond,
                       const Vector<SIMD<T1>, l> &if_true_val, const T2 &if_false_val) {
  for (int i = 0; i < l; ++i) r_val[i] = cond[i] ? if_true_val[i] : if_false_val;
}

// scalar, vector

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_MISC_RETURN_TYPE(T, l)
    conditional_assign(const Vector<SIMD<bool>, l> &cond, const T &if_true_val,
                       const Vector<SIMD<T>, l> &if_false_val) {
  Vector<SIMD<T>, l> r_val;
  for (int i = 0; i < l; ++i) r_val[i] = cond[i] ? if_true_val : if_false_val[i];
  return r_val;
}

template <typename T0, typename T1, typename T2, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_MISC_CONVERTIBLE_RETURN_VOID_TYPE(T0, T1, T2, l)
    conditional_assign(/* */ Vector<SIMD<T0>, l> &r_val, const Vector<SIMD<bool>, l> &cond, const T1 &if_true_val,
                       const Vector<SIMD<T2>, l> &if_false_val) {
  for (int i = 0; i < l; ++i) r_val[i] = cond[i] ? if_true_val : if_false_val[i];
}

// vector, vector

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_MISC_RETURN_TYPE(T, l)
    conditional_assign(const Vector<SIMD<bool>, l> &cond, const Vector<SIMD<T>, l> &if_true_val,
                       const Vector<SIMD<T>, l> &if_false_val) {
  Vector<SIMD<T>, l> r_val;
  for (int i = 0; i < l; ++i) r_val[i] = cond[i] ? if_true_val[i] : if_false_val[i];
  return r_val;
}

template <typename T0, typename T1, typename T2, int l>
KOKKOS_INLINE_FUNCTION static KOKKOSKERNELS_SIMD_MISC_CONVERTIBLE_RETURN_VOID_TYPE(T0, T1, T2, l)
    conditional_assign(/* */ Vector<SIMD<T0>, l> &r_val, const Vector<SIMD<bool>, l> &cond,
                       const Vector<SIMD<T1>, l> &if_true_val, const Vector<SIMD<T2>, l> &if_false_val) {
  for (int i = 0; i < l; ++i) r_val[i] = cond[i] ? if_true_val[i] : if_false_val[i];
}

template <typename T, int l, typename BinaryOp>
KOKKOS_INLINE_FUNCTION static T reduce(const Vector<SIMD<T>, l> &val, const BinaryOp &func) {
  T r_val = val[0];
  for (int i = 1; i < l; ++i) r_val = func(r_val, val[i]);
  return r_val;
}

template <typename T, int l, typename BinaryOp>
KOKKOS_INLINE_FUNCTION static T reduce(const Vector<SIMD<T>, l> &val, const BinaryOp &func, const T init) {
  T r_val = init;
  for (int i = 0; i < l; ++i) r_val = func(r_val, val[i]);
  return r_val;
}

template <int l>
KOKKOS_INLINE_FUNCTION static bool is_all_true(const Vector<SIMD<bool>, l> &cond) {
  return reduce(cond, [](const bool left, const bool right) -> bool { return (left && right); });
}

template <int l>
KOKKOS_INLINE_FUNCTION static bool is_any_true(const Vector<SIMD<bool>, l> &cond) {
  return reduce(cond, [](const bool left, const bool right) -> bool { return left || right; });
}

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static T min(const Vector<SIMD<T>, l> &val) {
  return reduce(val, [](const T left, const T right) -> T {
    const auto tmp = left < right;
    return tmp * left + !tmp * right;
  });
}

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static T max(const Vector<SIMD<T>, l> &val) {
  return reduce(val, [](const T left, const T right) -> T {
    const auto tmp = left > right;
    return tmp * left + !tmp * right;
  });
}

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static T sum(const Vector<SIMD<T>, l> &val) {
  return reduce(
      val, [](const T left, const T right) -> T { return left + right; }, T(0));
}

template <typename T, int l>
KOKKOS_INLINE_FUNCTION static T prod(const Vector<SIMD<T>, l> &val) {
  return reduce(
      val, [](const T left, const T right) -> T { return left * right; }, T(1));
}

#undef KOKKOSKERNELS_SIMD_MISC_RETURN_TYPE
#undef KOKKOSKERNELS_SIMD_MISC_CONVERTIBLE_RETURN_VOID_TYPE

}  // namespace KokkosBatched

#endif
