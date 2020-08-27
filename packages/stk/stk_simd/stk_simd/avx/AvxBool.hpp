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
// 

#ifndef STK_SIMD_AVX_BOOLS_H
#define STK_SIMD_AVX_BOOLS_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <immintrin.h>

namespace stk {
namespace simd {

class Bool {

 public:
    
  STK_MATH_FORCE_INLINE Bool() {}
    
  STK_MATH_FORCE_INLINE Bool(const bool x)
    : _data(x ? _mm256_cmp_pd(_mm256_setzero_pd(),_mm256_setzero_pd(),_CMP_EQ_OQ) : _mm256_setzero_pd())
    {
    }
    
  STK_MATH_FORCE_INLINE Bool(const __m256d& x)
    : _data(x)
    {
    }

  STK_MATH_FORCE_INLINE Bool(const Bool& x)
    : _data(x._data)
    {
    }
    
  STK_MATH_FORCE_INLINE Bool& operator= (const Bool& x) {
    _data = x._data;
    return *this;
  }
    
  STK_MATH_FORCE_INLINE double& operator[](int i) {return (reinterpret_cast<double*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const double& operator[](int i) const {return (reinterpret_cast<const double*>(&_data))[i];}

  __m256d _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};
  
} // namespace simd
} // namespace stk

#endif // SIMD_AVX_BOOLS_H__
