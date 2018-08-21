// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_AVX512_BOOLSF_H
#define STK_SIMD_AVX512_BOOLSF_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <immintrin.h>

namespace stk {
namespace simd {

class Boolf {

 public:
    
  STK_MATH_FORCE_INLINE Boolf() {}
    
  STK_MATH_FORCE_INLINE Boolf(const bool x) 
    : _data(x ? 
            _mm512_cmp_ps_mask(_mm512_setzero_ps(),_mm512_setzero_ps(),_CMP_EQ_OQ) :
            _mm512_cmp_ps_mask(_mm512_setzero_ps(),_mm512_setzero_ps(),_CMP_LT_OS) )
    {
    }
    
  STK_MATH_FORCE_INLINE Boolf(const __mmask16& x) 
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

  STK_MATH_FORCE_INLINE double operator[](int i) const {
    __m512 tmp = _mm512_mask_blend_ps(_data,Float(0.0)._data,Float(1.0)._data);
    return (reinterpret_cast<const double*>(&tmp))[i];
  }

  __mmask16 _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

const Boolf TRUE_VALf(true);
const Boolf FALSE_VALf(false);

} // namespace simd
} // namespace stk

#endif // SIMD_AVX512_BOOLSF_H__
