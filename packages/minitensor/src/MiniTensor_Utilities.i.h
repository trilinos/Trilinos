// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Utilities_i_h)
#define MiniTensor_Utilities_i_h

#include <cfloat>
#include <cmath>
#include <limits>
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace minitensor {

//
//
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
abs(T const & a)
{
  return a < T(0) ? -a : a;
}

//
//
//
template<typename T>
KOKKOS_INLINE_FUNCTION
void
swap(T & a, T & b)
{
  // Guard against the same memory location.
  if (&a == &b) return;

  auto const c = a;
  a = b;
  b = c;
}

//
//
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
max(T const & a, T const & b)
{
  return a > b ? a : b;
}

//
//
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
min(T const & a, T const & b)
{
  return a < b ? a : b;
}

//
// Sign function
//
template <typename T>
KOKKOS_INLINE_FUNCTION
int
sgn(T const & s)
{
  return (T(0) < s) - (s < T(0));
}

//
// Copysign function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
copysign(T const & a, T const & b)
{
  return b >= 0 ? T(std::abs(a)) : T(-std::abs(a));
}

//
// NaN function. Necessary to choose the proper underlying NaN
// for non-floating-point types.
// Assumption: non-floating-point types have a typedef that
// determines the underlying floating-point type.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
not_a_number()
{
  using S = typename Sacado::ScalarType<T>::type;
  return Kokkos::ArithTraits<S>::nan();
}

//
// Machine epsilon function. Necessary to choose the proper underlying
// machine epsilon for non-floating-point types.
// Assumption: non-floating-point types have a typedef that
// determines the underlying floating-point type.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
machine_epsilon()
{
  using S = typename Sacado::ScalarType<T>::type;
  return Kokkos::ArithTraits<S>::epsilon();
}

//
// Number of digits for integer types.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
Index
num_digits()
{
  return 0;
}

template<>
KOKKOS_INLINE_FUNCTION
Index
num_digits<Index>()
{
  return INDEX_SIZE;
}

template<>
KOKKOS_INLINE_FUNCTION
Index
num_digits<LongIndex>()
{
  return LONG_INDEX_SIZE;
}

//
// The circle constant
//
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
tau()
{
  using S = typename Sacado::ScalarType<T>::type;
  return static_cast<S>(
      6.283185307179586476925286766559005768394338798750211641949889185
  );
}

//
// Random number generation. Teuchos [-1,1]
//
template <typename T> typename Sacado::ScalarType<T>::type random() {
  using S = typename Sacado::ScalarType<T>::type;
  return Teuchos::ScalarTraits<S>().random();
}

//
// Uniform [0,1] random number generation.
//
template <typename T> typename Sacado::ScalarType<T>::type random_uniform() {
  using S = typename Sacado::ScalarType<T>::type;
  return static_cast<S>(0.5 * random<S>() + 0.5);
}

//
// Normal N(0,1) random number generation.
//
template <typename T> typename Sacado::ScalarType<T>::type random_normal() {
  using S = typename Sacado::ScalarType<T>::type;

  S const
  R = random_uniform<S>();

  S const
  Theta = tau<S>() * random_uniform<S>();
  return static_cast<S>(std::sqrt(-2.0 * std::log(R)) * cos(Theta));
}

//
// Fill in all levels of AD with specified constant.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
void
fill_AD(
    typename enable_if<is_same<T, typename ScalarType<T>::type>::value, T>::type & x,
    typename ScalarType<T>::type const c)
{
  x = c;
  return;
}

template<typename T>
KOKKOS_INLINE_FUNCTION
void
fill_AD(
    typename enable_if<!is_same<T, typename ScalarType<T>::type>::value, T>::type & x,
    typename ScalarType<T>::type const c)
{
  auto const
  order = x.size();

  // No AD info. Nothing to do.
  if (order == 0) return;

  using S = typename Sacado::ValueType<T>::type;

  for (auto i = 0; i < order; ++i) {
    fill_AD<S>(x.fastAccessDx(i), c);
  }

  return;
}

//
// Compute a non-negative integer power by binary manipulation.
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
integer_power(T const & X, Index const exponent)
{
  if (X == 0 || X == 1) return X;

  switch (exponent) {
    default:
      break;
    case 0:
      return 1;
    case 1:
      return X;
    case 2:
      return X * X;
    case 3:
      return X * X * X;
    case 4:
    {
      T const Y = X * X;
      return Y * Y;
    }
  }

  Index const
  rightmost_bit = 1;

  Index const
  number_digits = num_digits<Index>();

  Index const
  leftmost_bit = rightmost_bit << (number_digits - 1);

  Index
  t = 0;

  for (Index j = 0; j < number_digits; ++j) {

    if (((exponent << j) & leftmost_bit) != 0) {

      t = number_digits - j - 1;
      break;

    }

  }

  T
  P = X;

  Index
  i = 0;

  Index
  m = exponent;

  while ((m & rightmost_bit) == 0) {
    P = P * P;
    ++i;
    m = m >> 1;
  }

  T
  Y = P;

  for (Index j = i + 1; j <= t; ++j) {
    P = P * P;

    if (((exponent >> j) & rightmost_bit) != 0) {
      Y = Y * P;
    }
  }

  return Y;
}

//
// Integer nth root
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
integer_root(T const & x, Index const root)
{
  assert(root > 0);
  assert(x >= 0);

  if (root == 1 || x == 0 || x == 1) return x;

  T hi = 1;

  while (integer_power(hi, root) < x) hi *= 2;

  T lo = hi / 2;

  while (hi - lo > 1) {
    T mid = (lo + hi) / 2;
    T t = integer_power(mid, root);
    if (t < x) lo = mid;
    else if (x < t) hi = mid;
    else return mid;
  }

  if (integer_power(hi, root) == x) return hi;

  return lo;
}

//
// Utility for Kronecker delta in 2D
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j)
{
  assert(0 <= i && i < 2);
  assert(0 <= j && j < 2);

  if (i == j) return T(1);

  return T(0);
}

//
// Utility for Kronecker delta in 3D
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j, Index const k)
{
  assert(0 <= i && i < 3);
  assert(0 <= j && j < 3);
  assert(0 <= k && k < 3);

  if (i == j && j == k) return T(1);

  return T(0);
}

//
// Utility for Kronecker delta in 4D
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j, Index const k, Index const l)
{
  assert(0 <= i && i < 4);
  assert(0 <= j && j < 4);
  assert(0 <= k && k < 4);
  assert(0 <= l && l < 4);

  if (i == j && j == k && k == l) return T(1);

  return T(0);
}

//
// Utility for Levi-Civita/permutation/alternating symbol in 2D
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j)
{
  assert(0 <= i && i < 2);
  assert(0 <= j && j < 2);

  if (i == 0 && j == 1) return T(1);

  if (i == 1 && j == 0) return T(-1);

  return T(0);
}

//
// Utility for Levi-Civita/permutation/alternating symbol in 3D
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j, Index const k)
{
  assert(0 <= i && i < 3);
  assert(0 <= j && j < 3);
  assert(0 <= k && k < 3);

  if (i == 0 && j == 1 && k == 2) return T(1);
  if (i == 1 && j == 2 && k == 0) return T(1);
  if (i == 2 && j == 0 && k == 1) return T(1);

  if (i == 2 && j == 1 && k == 0) return T(-1);
  if (i == 0 && j == 2 && k == 1) return T(-1);
  if (i == 1 && j == 0 && k == 2) return T(-1);

  return T(0);
}

//
// Utility for Levi-Civita/permutation/alternating symbol in 4D
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j, Index const k, Index const l)
{
  assert(0 <= i && i < 4);
  assert(0 <= j && j < 4);
  assert(0 <= k && k < 4);
  assert(0 <= l && l < 4);

  if (i == 0 && j == 1 && k == 2 && l == 3) return T(1);
  if (i == 1 && j == 2 && k == 3 && l == 0) return T(1);
  if (i == 2 && j == 3 && k == 0 && l == 1) return T(1);
  if (i == 3 && j == 0 && k == 1 && l == 2) return T(1);

  if (i == 3 && j == 2 && k == 1 && l == 0) return T(-1);
  if (i == 0 && j == 3 && k == 2 && l == 1) return T(-1);
  if (i == 1 && j == 0 && k == 3 && l == 2) return T(-1);
  if (i == 2 && j == 1 && k == 0 && l == 3) return T(-1);

  return T(0);
}

} // namespace minitensor

#endif // MiniTensor_Utilities_i_h
