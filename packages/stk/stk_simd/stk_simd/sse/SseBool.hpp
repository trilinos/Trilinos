// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_SSE_BOOLS_H
#define STK_SIMD_SSE_BOOLS_H

#include <emmintrin.h>

namespace stk {

namespace simd {

class Bool {

 public:
    
  inline Bool() : _data(_mm_setzero_pd()) {}
    
  inline Bool(const bool x) 
    : _data(x ? _mm_cmpeq_pd(_mm_setzero_pd(),_mm_setzero_pd()) : _mm_setzero_pd())
    {
    }
    
  inline Bool(const __m128d& x) 
    : _data(x) 
    {
    }

  inline Bool(const Bool& x) 
    : _data(x._data) 
    {
    }
    
  inline Bool& operator= (const Bool& x) {
    _data = x._data;
    return *this;
  }

  inline double& operator[](int i) {return (reinterpret_cast<double*>(&_data))[i];}
  inline const double& operator[](int i) const {return (reinterpret_cast<const double*>(&_data))[i];}
     
  __m128d _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd

} // namespace stk

#endif // SIMD_SSE_BOOLS_H__
