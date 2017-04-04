// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_AVX512_BOOLSF_H
#define STK_SIMD_AVX512_BOOLSF_H

#include <immintrin.h>

namespace stk {
namespace simd {

class Boolf {

 public:
    
  inline Boolf() : _data() {}
    
  inline Boolf(const bool x) 
    : _data(x ? 
            _mm512_cmp_ps_mask(_mm512_setzero_ps(),_mm512_setzero_ps(),_CMP_EQ_OQ) :
            _mm512_cmp_ps_mask(_mm512_setzero_ps(),_mm512_setzero_ps(),_CMP_LT_OS) )
    {
    }
    
  inline Boolf(const __mmask16& x) 
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

  inline double operator[](int i) const {
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
