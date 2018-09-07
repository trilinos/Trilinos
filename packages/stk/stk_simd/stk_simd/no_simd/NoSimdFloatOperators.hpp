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

#include "./NoSimdSizes.hpp"

namespace stk {
namespace simd {

namespace hidden {
static const simd::Boolf TRUE_VALf(true);
}

STK_MATH_FORCE_INLINE simd::Float operator+ (const simd::Float& a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] + b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator- (const simd::Float& a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] - b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator* (const simd::Float& a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] * b[i];
  return tmp;
}
   
STK_MATH_FORCE_INLINE simd::Float operator/ (const simd::Float& a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] / b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator+ (const float a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a + b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator+ (const simd::Float& a, const float b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] + b;
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator- (const float a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a - b[i];
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Float operator- (const simd::Float& a, const float b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] - b;
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator* (const float a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a * b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator* (const simd::Float& a, const float b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] * b;
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float operator/ (const float a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a / b[i];
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Float operator/ (const simd::Float& a, const float b) {
  simd::Float tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] / b;
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator< (const simd::Float& a, const simd::Float& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] < b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator<=(const simd::Float& a, const simd::Float& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] <= b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator> (const simd::Float& a, const simd::Float& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] > b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator>=(const simd::Float& a, const simd::Float& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] >= b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator==(const simd::Float& a, const simd::Float& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] == b[i];
  return tmp;
}
  
STK_MATH_FORCE_INLINE simd::Boolf operator!=(const simd::Float& a, const simd::Float& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] != b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator&& (const simd::Boolf& a, const simd::Boolf& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] && b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator|| (const simd::Boolf& a, const simd::Boolf& b) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = a[i] || b[i];
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Boolf operator! (const simd::Boolf& a) {
  simd::Boolf tmp;
  for (int i=0; i < nfloats; ++i) tmp[i] = ! a[i];
  return tmp;
}

} // namespace simd
} // namespace stk
