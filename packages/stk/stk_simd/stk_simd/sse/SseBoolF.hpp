// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_SSE_BOOLSF_H
#define STK_SIMD_SSE_BOOLSF_H

#include <emmintrin.h>

namespace stk {
namespace simd {

class Boolf {

 public:
    
  inline Boolf() : _data(_mm_setzero_ps()) {}
    
  inline Boolf(const bool x) 
    : _data(x ? _mm_cmpeq_ps(_mm_setzero_ps(),_mm_setzero_ps()) : _mm_setzero_ps())
    {
    }
    
  inline Boolf(const __m128& x) 
    : _data(x) 
    {
    }

  inline Boolf(const Boolf& x) 
    : _data(x._data) 
    {
    }
    
  inline Boolf& operator= (const Boolf& x) {
    _data = x._data;
    return *this;
  }

  inline float& operator[](int i) {return (reinterpret_cast<float*>(&_data))[i];}
  inline const float& operator[](int i) const {return (reinterpret_cast<const float*>(&_data))[i];}
     
  __m128 _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd
} // namespace stk

#endif // SIMD_SSE_BOOLSF_H__
