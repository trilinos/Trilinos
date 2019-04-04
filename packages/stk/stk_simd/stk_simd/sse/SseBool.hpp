// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_SSE_BOOLS_H
#define STK_SIMD_SSE_BOOLS_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <emmintrin.h>

namespace stk {
namespace simd {

class Bool {

 public:
    
  STK_MATH_FORCE_INLINE Bool() {}
    
  STK_MATH_FORCE_INLINE Bool(const bool x) 
    : _data(x ? _mm_cmpeq_pd(_mm_setzero_pd(),_mm_setzero_pd()) : _mm_setzero_pd())
    {
    }
    
  STK_MATH_FORCE_INLINE Bool(const __m128d& x) 
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
     
  __m128d _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd

} // namespace stk

#endif // SIMD_SSE_BOOLS_H__
