// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_DEFAULT_BOOLS_H
#define STK_SIMD_DEFAULT_BOOLS_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <array>
#include "./NoSimdSizes.hpp"

namespace stk {
namespace simd {

class Bool {

 public:
    
  STK_MATH_FORCE_INLINE Bool() {}
    
  STK_MATH_FORCE_INLINE Bool(const bool x) {
    for (int i=0; i < ndoubles; ++i) _data[i] = x;
  }
    
  STK_MATH_FORCE_INLINE Bool(const Bool& x) {
    for (int i=0; i < ndoubles; ++i) _data[i] = x[i];
  }
    
  STK_MATH_FORCE_INLINE Bool& operator= (const Bool& x) {
    for (int i=0; i < ndoubles; ++i) _data[i] = x[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE bool& operator[](int i) {return _data[i];}
  STK_MATH_FORCE_INLINE const bool& operator[](int i) const {return _data[i];}
     
  bool _data[ndoubles];
};

} // namespace simd
} // namespace stk

#endif // STK_SIMD_DEFAULT_BOOLS_H__
