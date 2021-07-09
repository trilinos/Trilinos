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

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <Kokkos_Macros.hpp>

#define STK_MATH_FORCE_INLINE KOKKOS_FORCEINLINE_FUNCTION

namespace stk {
namespace math {

namespace hidden {
static STK_MATH_FORCE_INLINE double Cbrt(double x) { return cbrt(x); }
static STK_MATH_FORCE_INLINE float CbrtF(float x) { return cbrtf(x); }
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T fmadd(T a, T b, T c) {
  return a*b+c;
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T sqrt(T x) {
  return std::sqrt(x);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T cbrt(T x) {
  return hidden::Cbrt(x);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T log(T x) {
  return std::log(x);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T log10(T x) {
  return std::log10(x);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T exp(T x) {
  return std::exp(x);
}
   
STK_MATH_FORCE_INLINE double pow(const double x, const double y) {
  return std::pow(x, y);
}

STK_MATH_FORCE_INLINE float pow(const float x, const float y) {
  return std::pow(x, y);
}

STK_MATH_FORCE_INLINE double pow(const double x, int y) {
  return std::pow(x, y);
}

STK_MATH_FORCE_INLINE float pow(const float x, int y) {
  return std::pow(x, y);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T sin(T a) {
  return std::sin(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T cos(T a) {
  return std::cos(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T tan(T a) {
  return std::tan(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T sinh(T a) {
  return std::sinh(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T cosh(T a) {
  return std::cosh(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T tanh(T a) {
  return std::tanh(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T asin(T a) {
  return std::asin(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T acos(T a) {
  return std::acos(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T atan(T a) {
  return std::atan(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T atan2(T a, T b) {
  return std::atan2(a, b);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T asinh(T a) {
  return std::asinh(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T acosh(T a) {
  return std::acosh(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T atanh(T a) {
  return std::atanh(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T erf(T a) {
  return std::erf(a);
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T multiplysign(T x, T y) { // return x times sign of y
  return x * std::copysign(static_cast<T>(1.0), y);
}

//template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
//STK_MATH_FORCE_INLINE T copysign(T x, T y) { // return abs(x) times sign of y
//  return std::copysign(x, y);
//}

STK_MATH_FORCE_INLINE double copysign(double x, double y) { // return abs(x) times sign of y
  return std::copysign(x, y);
}

STK_MATH_FORCE_INLINE float copysign(float x, float y) { // return abs(x) times sign of y
  return std::copysign(x, y);
}

STK_MATH_FORCE_INLINE double abs(const double x) {
  return std::abs(x);
}

STK_MATH_FORCE_INLINE float abs(const float x) {
  return std::abs(x);
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
  return std::isnan(a);
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

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
STK_MATH_FORCE_INLINE T if_then_else_zero(const bool b, T v) {
  return b ? v : static_cast<T>(0);
}

} // namespace math
} // namespace stk

#endif
