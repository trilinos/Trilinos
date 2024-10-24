// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef STK_MATH_H
#define STK_MATH_H

#include <cstdlib>
#include <cassert>
#include <type_traits>
#include <Kokkos_Macros.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

#define STK_MATH_FORCE_INLINE KOKKOS_FORCEINLINE_FUNCTION

namespace stk {
namespace math {

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T fmadd(T a, T b, T c) {
  return a*b+c;
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T sqrt(T x)
{
  return Kokkos::sqrt(x);
}

inline long double sqrt(long double x) {
  return Kokkos::sqrtl(x);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T cbrt(T x) {
  return Kokkos::cbrt(x);
}

inline long double cbrt(long double x) {
  return Kokkos::cbrtl(x);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T log(T x) {
  return Kokkos::log(x);
}

inline long double log(long double x) {
  return Kokkos::logl(x);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T log10(T x) {
  return Kokkos::log10(x);
}

inline long double log10(long double x) {
  return Kokkos::log10l(x);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T exp(T x) {
  return Kokkos::exp(x);
}

inline long double exp(long double x) {
  return Kokkos::expl(x);
}
   
STK_MATH_FORCE_INLINE double pow(const double x, const double y) {
  return Kokkos::pow(x, y);
}

STK_MATH_FORCE_INLINE float pow(const float x, const float y) {
  return Kokkos::pow(x, y);
}

STK_MATH_FORCE_INLINE double pow(const double x, int y) {
  return Kokkos::pow(x, y);
}

STK_MATH_FORCE_INLINE float pow(const float x, int y) {
  return Kokkos::pow(x, y);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T sin(T a) {
  return Kokkos::sin(a);
}

inline long double sin(long double a) {
  return Kokkos::sinl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T cos(T a) {
  return Kokkos::cos(a);
}

inline long double cos(long double a) {
  return Kokkos::cosl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T tan(T a) {
  return Kokkos::tan(a);
}

inline long double tan(long double a) {
  return Kokkos::tanl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T sinh(T a) {
  return Kokkos::sinh(a);
}

inline long double sinh(long double a) {
  return Kokkos::sinhl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T cosh(T a) {
  return Kokkos::cosh(a);
}

inline long double cosh(long double a) {
  return Kokkos::coshl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T tanh(T a) {
  return Kokkos::tanh(a);
}

inline long double tanh(long double a) {
  return Kokkos::tanhl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T asin(T a) {
  return Kokkos::asin(a);
}

inline long double asin(long double a) {
  return Kokkos::asinl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T acos(T a) {
  return Kokkos::acos(a);
}

inline long double acos(long double a) {
  return Kokkos::acosl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T atan(T a) {
  return Kokkos::atan(a);
}

inline long double atan(long double a) {
  return Kokkos::atanl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T atan2(T a, T b) {
  return Kokkos::atan2(a,b);
}

inline long double atan2(long double a, long double b) {
  return Kokkos::atan2l(a,b);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T asinh(T a) {
  return Kokkos::asinh(a);
}

inline long double asinh(long double a) {
  return Kokkos::asinhl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T acosh(T a) {
  return Kokkos::acosh(a);
}

inline long double acosh(long double a) {
  return Kokkos::acoshl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T atanh(T a) {
  return Kokkos::atanh(a);
}

inline long double atanh(long double a) {
  return Kokkos::atanhl(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T erf(T a) {
  return Kokkos::erf(a);
}

inline long double erf(long double a) {
  return Kokkos::erf(a);
}

template <typename T,
    typename std::enable_if<std::is_floating_point<T>::value && !std::is_same<T, long double>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T copysign(T x, T y) {
  return Kokkos::copysign(x,y);
}

inline long double copysign(long double x, long double y) {
  return Kokkos::copysignl(x,y);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T multiplysign(T x, T y) { // return x times sign of y
  return x * copysign(static_cast<T>(1.0), y);
}

STK_MATH_FORCE_INLINE double abs(const double x) {
  return Kokkos::abs(x);
}

STK_MATH_FORCE_INLINE float abs(const float x) {
  return Kokkos::abs(x);
}

STK_MATH_FORCE_INLINE double max(const double x, const double y) {
  return x > y ? x : y;
}

STK_MATH_FORCE_INLINE float max(const float x, const float y) {
  return x > y ? x : y;
}

STK_MATH_FORCE_INLINE double min(const double& x, const double& y) {
  return x < y ? x : y;
}

STK_MATH_FORCE_INLINE float min(const float& x, const float& y) {
  return x < y ? x : y;
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE bool isnan(T a) {
  return Kokkos::isnan(a);
}

//template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
//STK_MATH_FORCE_INLINE T if_then_else(const bool b, T v1, T v2) {
//  return b ? v1 : v2;
//}

STK_MATH_FORCE_INLINE float if_then_else(const bool b, float v1, float v2) {
  return b ? v1 : v2;
}

STK_MATH_FORCE_INLINE double if_then_else(const bool b, double v1, double v2) {
  return b ? v1 : v2;
}

STK_MATH_FORCE_INLINE long double if_then_else(const bool b, long double v1, long double v2) {
  return b ? v1 : v2;
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T if_then_else_zero(const bool b, T v) {
  return b ? v : static_cast<T>(0);
}

} // namespace math
} // namespace stk

#endif
