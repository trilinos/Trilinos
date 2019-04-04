// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <stk_math/StkMath.hpp>
#include "NoSimdSizes.hpp"

namespace stk {
namespace math {

STK_MATH_FORCE_INLINE simd::Double fmadd(const simd::Double& a, const simd::Double& b, const simd::Double& c) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::fmadd(a[i], b[i], c[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double sqrt(const simd::Double& x) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::sqrt(x[i]);
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Double cbrt(const simd::Double& x) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::cbrt(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double log(const simd::Double& x) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::log(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double log10(const simd::Double& x) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::log10(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double exp(const simd::Double& x) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::exp(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const simd::Double& y) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::pow(x[i], y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const double y) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::pow(x[i], y);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const int y) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::pow(x[i], y);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double sin(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::sin(a[i]);
  return tmp;
}
 
STK_MATH_FORCE_INLINE simd::Double cos(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::cos(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double tan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::tan(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double sinh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::sinh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double cosh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::cosh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double tanh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::tanh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double asin(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::asin(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double acos(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::acos(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double atan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::atan(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double atan2(const simd::Double& a, const simd::Double& b) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::atan2(a[i],b[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double asinh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::asinh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double acosh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::acosh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double atanh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::atanh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double erf(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = std::erf(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double multiplysign(const simd::Double& x, const simd::Double& y) { // return x times sign of y
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::multiplysign(x[i], y[i]);
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Double copysign(const simd::Double& x, const simd::Double& y) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::copysign(x[i], y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double abs(const simd::Double& x) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::abs(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double min(const simd::Double& x, const simd::Double& y) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::min(x[i],y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double max(const simd::Double& x, const simd::Double& y) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::max(x[i],y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Bool isnan(const simd::Double& a) {
  simd::Bool tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::isnan(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double if_then_else(const simd::Bool& b, const simd::Double& v1, const simd::Double& v2) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::if_then_else(b[i], v1[i], v2[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double if_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) tmp[i] = stk::math::if_then_else_zero(b[i], v[i]);
  return tmp;
}

} // namespace math
} // namespace stk

