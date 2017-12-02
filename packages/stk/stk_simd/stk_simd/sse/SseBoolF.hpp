// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_SSE_BOOLSF_H
#define STK_SIMD_SSE_BOOLSF_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <emmintrin.h>

namespace stk {
namespace simd {

class Boolf {

 public:
    
  STK_MATH_FORCE_INLINE Boolf() {}
    
  STK_MATH_FORCE_INLINE Boolf(const bool x) 
    : _data(x ? _mm_cmpeq_ps(_mm_setzero_ps(),_mm_setzero_ps()) : _mm_setzero_ps())
    {
    }
    
  STK_MATH_FORCE_INLINE Boolf(const __m128& x) 
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

  STK_MATH_FORCE_INLINE float& operator[](int i) {return (reinterpret_cast<float*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const float& operator[](int i) const {return (reinterpret_cast<const float*>(&_data))[i];}
     
  __m128 _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd
} // namespace stk

#endif // SIMD_SSE_BOOLSF_H__
