// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Utilities_h)
#define MiniTensor_Utilities_h

#include <type_traits>

#include "MiniTensor_config.h"
#include "Sacado.hpp"

namespace minitensor {

//
// abs function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
abs(T const & a);

//
//swap function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
void
swap(T & a, T & b);

//
// max function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
max(T const & a, T const & b);

//
// max function
//
template<typename T>
KOKKOS_INLINE_FUNCTION
T
min(T const & a,T const & b);

///
/// Sign function
///
template<typename T>
KOKKOS_INLINE_FUNCTION
int
sgn(T const & s);

///
/// Copysign function
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
copysign(T const & a, T const & b);

///
/// NaN function. Necessary to choose the proper underlying NaN
/// for non-floating-point types.
/// Assumption: non-floating-point types have a typedef that
/// determines the underlying floating-point type.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
not_a_number();

///
/// Machine epsilon function. Necessary to choose the proper underlying
/// machine epsilon for non-floating-point types.
/// Assumption: non-floating-point types have a typedef that
/// determines the underlying floating-point type.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
machine_epsilon();

///
/// Number of digits for integer types.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
Index
num_digits();

template<>
KOKKOS_INLINE_FUNCTION
Index
num_digits<Index>();

template<>
KOKKOS_INLINE_FUNCTION
Index
num_digits<LongIndex>();

///
/// The circle constant
///
template<typename T>
KOKKOS_INLINE_FUNCTION
typename Sacado::ScalarType<T>::type
tau();

///
/// Random number generation. Uniform distribution U(-1,1)
/// which is the Teuchos default (!).
///
template <typename T> typename Sacado::ScalarType<T>::type random();

///
/// Random number generation. Uniform distribution U(0,1).
///
template <typename T> typename Sacado::ScalarType<T>::type random_uniform();

///
/// Random number generation. Normal distribution N(0,1).
///
template <typename T> typename Sacado::ScalarType<T>::type random_normal();

///
/// Fill all levels of AD to specified constant.
///
using std::enable_if;
using std::is_same;
using Sacado::ScalarType;

template<typename T>
KOKKOS_INLINE_FUNCTION
void
fill_AD(
    typename enable_if<is_same<T, typename ScalarType<T>::type>::value, T>::type & x,
    typename ScalarType<T>::type const c);

template<typename T>
KOKKOS_INLINE_FUNCTION
void
fill_AD(
    typename enable_if<!is_same<T, typename ScalarType<T>::type>::value, T>::type & x,
    typename ScalarType<T>::type const c);

///
/// Compute a non-negative integer power by binary manipulation.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
integer_power(T const & X, Index const exponent);

///
/// Compute a non-negative integer root by bisection.
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
integer_root(T const & x, Index const root);

///
/// Utility for Kronecker delta in 2D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j);

///
/// Utility for Kronecker delta in 3D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j, Index const k);

///
/// Utility for Kronecker delta in 4D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
kronecker_delta(Index const i, Index const j, Index const k, Index const l);

///
/// Utility for Levi-Civita/permutation/alternating symbol in 2D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j);

///
/// Utility for Levi-Civita/permutation/alternating symbol in 3D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j, Index const k);

///
/// Utility for Levi-Civita/permutation/alternating symbol in 4D
///
template<typename T>
KOKKOS_INLINE_FUNCTION
T
levi_civita(Index const i, Index const j, Index const k, Index const l);

} // namespace minitensor

#include "MiniTensor_Utilities.i.h"

#endif // MiniTensor_Utilities_h
