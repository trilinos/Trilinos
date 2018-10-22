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

namespace stk {
namespace simd {

namespace hidden {
static const simd::Boolf TRUE_VALf(true);
}

inline simd::Float operator+ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_add_ps(a._data,b._data));
}

inline simd::Float operator- (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_sub_ps(a._data,b._data));
}

inline simd::Float operator* (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_mul_ps(a._data,b._data));
}
   
inline simd::Float operator/ (const simd::Float& a, const simd::Float& b) {
  return simd::Float(_mm_div_ps(a._data,b._data));
}

inline simd::Float operator+ (const float a, const simd::Float& b) {
  return simd::Float(_mm_add_ps(_mm_set1_ps(a),b._data));
}

inline simd::Float operator+ (const simd::Float& a, const float b) {
  return simd::Float(_mm_add_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Float operator- (const float a, const simd::Float& b) {
  return simd::Float(_mm_sub_ps(_mm_set1_ps(a),b._data));
}
  
inline simd::Float operator- (const simd::Float& a, const float b) {
  return simd::Float(_mm_sub_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Float operator* (const float a, const simd::Float& b) {
  return simd::Float(_mm_mul_ps(_mm_set1_ps(a),b._data));
}

inline simd::Float operator* (const simd::Float& a, const float b) {
  return simd::Float(_mm_mul_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Float operator/ (const float a, const simd::Float& b) {
  return simd::Float(_mm_div_ps(_mm_set1_ps(a),b._data));
}
  
inline simd::Float operator/ (const simd::Float& a, const float b) {
  return simd::Float(_mm_div_ps(a._data,_mm_set1_ps(b)));
}

inline simd::Boolf operator< (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmplt_ps(x._data,y._data));
}

inline simd::Boolf operator<=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmple_ps(x._data,y._data));
}

inline simd::Boolf operator> (const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmplt_ps(y._data,x._data));
}

inline simd::Boolf operator>=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmple_ps(y._data,x._data));
}

inline simd::Boolf operator==(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmpeq_ps(x._data,y._data));
}
  
inline simd::Boolf operator!=(const simd::Float& x, const simd::Float& y) {
  return simd::Boolf(_mm_cmpneq_ps(x._data,y._data));
}

inline simd::Boolf operator&& (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf(_mm_and_ps(x._data,y._data));
}

inline simd::Boolf operator|| (const simd::Boolf& x, const simd::Boolf& y) {
  return simd::Boolf(_mm_or_ps(x._data,y._data));
}

inline simd::Boolf operator! (const simd::Boolf& x) {
  return simd::Boolf( _mm_xor_ps(hidden::TRUE_VALf._data,x._data) );
}

} // namespace simd
} // namespace stk
