// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_AVX512_BOOLS_H
#define STK_SIMD_AVX512_BOOLS_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <immintrin.h>

namespace stk {
namespace simd {

class Bool {

 public:
    
  STK_MATH_FORCE_INLINE Bool() {}
    
  STK_MATH_FORCE_INLINE  Bool(const bool x)
    : _data(x ? 
            _mm512_cmp_pd_mask(_mm512_setzero_pd(),_mm512_setzero_pd(),_CMP_EQ_OQ) : 
            _mm512_cmp_pd_mask(_mm512_setzero_pd(),_mm512_setzero_pd(),_CMP_LT_OQ) )
  {
  }
    
  STK_MATH_FORCE_INLINE  Bool(const __mmask8& x)
    : _data(x)
  {
  }
    
  STK_MATH_FORCE_INLINE Bool(const Bool& x)
    : _data(x._data)
  {
  }
    
  STK_MATH_FORCE_INLINE Bool& operator= (const Bool x) {
    _data = x._data;
    return *this;
  }
    
  // // These are for debugging only please
  // long long operator[](int i) {
  //   __m512i tmp = _mm512_broadcastmb_epi64(_data);
  //   return (reinterpret_cast<long long*>(&tmp))[i];
  // }

  STK_MATH_FORCE_INLINE double operator[](int i) const {
    __m512d tmp = _mm512_mask_blend_pd(_data,Double(0.0)._data,Double(1.0)._data);
    return (reinterpret_cast<const double*>(&tmp))[i];
  }

  __mmask8 _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};
  
} // namespace simd

} // namespace stk

#endif // SIMD_AVX512_BOOLS_H__
