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

STK_MATH_FORCE_INLINE simd::Float load_aligned(const float* x) {
  return simd::Float(_mm256_load_ps(x));
}

STK_MATH_FORCE_INLINE simd::Float load(const float* x) {
  return simd::Float(_mm256_loadu_ps(x));
}
    
STK_MATH_FORCE_INLINE simd::Float load(const float* x, const int offset) {
  return simd::Float(_mm256_setr_ps(x[0],       x[offset],  x[2*offset],x[3*offset],
                                     x[4*offset],x[5*offset],x[6*offset],x[7*offset]));
}

STK_MATH_FORCE_INLINE void store_aligned(float* x, const simd::Float& z) {
  _mm256_store_ps(x,z._data);
}

STK_MATH_FORCE_INLINE void store(float* x, const simd::Float& z) {
  _mm256_storeu_ps(x,z._data);
}
  
STK_MATH_FORCE_INLINE void store(float* x, const simd::Float& z, const int offset) {
  for (int i=0; i < nfloats; ++i) {
    x[offset*i] = z[i];
  }
}

} // namespace simd
} // namespace stk
