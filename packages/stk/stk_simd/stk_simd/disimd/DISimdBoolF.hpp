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

#ifndef STK_SIMD_DISIMD_DISIMDBOOLF_HPP
#define STK_SIMD_DISIMD_DISIMDBOOLF_HPP

namespace stk {
namespace simd {

class Boolf {

 public:
    
  STK_MATH_FORCE_INLINE Boolf() {}
    
  STK_MATH_FORCE_INLINE Boolf(bool x) 
    : _data(x)
    {
    }
    
  STK_MATH_FORCE_INLINE Boolf(const SIMD_NAMESPACE::simd_mask<float, SIMD_NAMESPACE::simd_abi::native>& x)
    : _data(x) 
    {
    }

  STK_MATH_FORCE_INLINE Boolf(const Boolf& x) 
    : _data(x._data) 
    {
    }
    
  STK_MATH_FORCE_INLINE Boolf& operator= (const Boolf& x) {
    _data = x._data;
    return *this;
  }

#if defined(__AVX512F__) && !defined(__CUDACC__) && !defined(STK_SIMD_NONE)
  STK_MATH_FORCE_INLINE float operator[](int i) const {
    __m512 tmp = _mm512_mask_blend_ps(_data.get(), _mm512_set1_ps(0.0), _mm512_set1_ps(1.0));
    return (reinterpret_cast<const float*>(&tmp))[i];
  }
#elif defined(__CUDACC__) || defined(STK_SIMD_NONE)
  STK_MATH_FORCE_INLINE float operator[](int i) const {
    return _data.get() ? 1.0f : 0.0f;
  }
#else
  STK_MATH_FORCE_INLINE float& operator[](int i) {return (reinterpret_cast<float*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const float& operator[](int i) const {return (reinterpret_cast<const float*>(&_data))[i];}
#endif
     
  SIMD_NAMESPACE::simd_mask<float, SIMD_NAMESPACE::simd_abi::native> _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

const Boolf TRUE_VALf(true);
const Boolf FALSE_VALf(false);

} // namespace simd
} // namespace stk

#endif // STK_SIMD_DISIMD_DISIMDBOOLF_HPP

