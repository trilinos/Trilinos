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
#include "./NoSimdSizes.hpp"

namespace stk {
namespace math {

STK_MATH_FORCE_INLINE simd::Float fmadd(const simd::Float& a, const simd::Float& b, const simd::Float& c) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::fmadd(a[i], b[i], c[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float sqrt(const simd::Float& x) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::sqrt(x[i]);
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Float cbrt(const simd::Float& x) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::cbrt(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float log(const simd::Float& x) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::log(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float log10(const simd::Float& x) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::log10(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float exp(const simd::Float& x) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::exp(x[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, const simd::Float& y) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::pow(x[i], y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, float y) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::pow(x[i], y);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, int y) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::pow(x[i], y);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float sin(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::sin(a[i]);
  return tmp;
}
 
STK_MATH_FORCE_INLINE simd::Float cos(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::cos(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float tan(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::tan(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float sinh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::sinh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float cosh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::cosh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float tanh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::tanh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float asin(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::asin(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float acos(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::acos(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float atan(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::atan(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float atan2(const simd::Float& a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::atan2(a[i],b[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float asinh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::asinh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float acosh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::acosh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float atanh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::atanh(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float erf(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = std::erf(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float multiplysign(const simd::Float& x, const simd::Float& y) { // return x times sign of y
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::multiplysign(x[i], y[i]);
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Float copysign(const simd::Float& x, const simd::Float& y) { // return abs(x) times sign of y
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::copysign(x[i], y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float abs(const simd::Float& x) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::abs(x[i]);
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Float min(const simd::Float& x, const simd::Float& y) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::min(x[i],y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float max(const simd::Float& x, const simd::Float& y) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::max(x[i],y[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf isnan(const simd::Float& a) {
  simd::Boolf tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::isnan(a[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float if_then_else(const simd::Boolf& b, const simd::Float& v1, const simd::Float& v2) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::if_then_else(b[i], v1[i], v2[i]);
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float if_then_else_zero(const simd::Boolf& b, const simd::Float& v) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) tmp[i] = stk::math::if_then_else_zero(b[i], v[i]);
  return tmp;
}

} // namespace math
} // namespace stk

