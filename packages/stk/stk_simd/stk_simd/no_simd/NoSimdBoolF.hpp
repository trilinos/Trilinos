// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_DEFAULT_BOOLSF_H
#define STK_SIMD_DEFAULT_BOOLSF_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <array>
#include "./NoSimdSizes.hpp"

namespace stk {
namespace simd {

class Boolf {

 public:
    
  STK_MATH_FORCE_INLINE Boolf() {}
    
  STK_MATH_FORCE_INLINE Boolf(const bool x) {
    for (int i=0; i < nfloats; ++i) _data[i] = x;   
  }

  STK_MATH_FORCE_INLINE Boolf(const Boolf& x) {
    for (int i=0; i < nfloats; ++i) _data[i] = x[i];
  }
  
  STK_MATH_FORCE_INLINE Boolf& operator= (const Boolf& x) {
    for (int i=0; i < nfloats; ++i) _data[i] = x[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE bool& operator[](int i) {return _data[i];}
  STK_MATH_FORCE_INLINE const bool& operator[](int i) const {return _data[i];}
     
  bool _data[nfloats];
};

} // namespace simd
} // namespace stk

#endif // STK_SIMD_DEFAULT_BOOLSF_H__
