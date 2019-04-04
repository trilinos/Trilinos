// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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

// IWYU pragma: private, include <stk_simd/Simd.hpp>

namespace stk {
namespace math {

namespace hidden {
static inline float CBRT(const float x) { return cbrtf(x); }
static const simd::Float SIGN_MASKf(-0);
}

inline simd::Float fmadd(const simd::Float& a, const simd::Float& b, const simd::Float& c) {
  return simd::Float(_mm_add_ps(_mm_mul_ps(a._data,b._data), c._data));
}

inline simd::Float sqrt(const simd::Float& x) {
  return simd::Float(_mm_sqrt_ps(x._data));
}
  
inline simd::Float cbrt(const simd::Float& x) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm_cbrt_ps(x._data));
#else
  simd::Float out;
  for (int i=0; i < simd::nfloats; ++i) {
    out[i] = hidden::CBRT(x[i]);
  }
  return out;
#endif
}

inline simd::Float log(const simd::Float& x) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm_log_ps(x._data));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::log(x[n]);
  }
  return tmp;
#endif
}

inline simd::Float log10(const simd::Float& x) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::log10(x[n]);
  }
  return tmp;
}

inline simd::Float exp(const simd::Float& x) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm_exp_ps(x._data));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::exp(x[n]);
  }
  return tmp;
#endif
}

inline simd::Float pow(const simd::Float& x, const simd::Float& y) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y[n]);
  }
  return tmp;
}

inline simd::Float pow(const simd::Float& x, float y) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
}

inline simd::Float pow(const simd::Float& x, int y) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
}

inline simd::Float sin(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::sin(a[i]);
  }
  return tmp;
}
 
inline simd::Float cos(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::cos(a[i]);
  }
  return tmp;
}

inline simd::Float tan(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::tan(a[i]);
  }
  return tmp;
}

inline simd::Float sinh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::sinh(a[i]);
  }
  return tmp;
}

inline simd::Float cosh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::cosh(a[i]);
  }
  return tmp;
}

inline simd::Float tanh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::tanh(a[i]);
  }
  return tmp;
}

inline simd::Float asin(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::asin(a[i]);
  }
  return tmp;
}

inline simd::Float acos(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::acos(a[i]);
  }
  return tmp;
}

inline simd::Float atan(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::atan(a[i]);
  }
  return tmp;
}

inline simd::Float atan2(const simd::Float& a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::atan2(a[i],b[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float asinh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::asinh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float acosh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::acosh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float atanh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::atanh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float erf(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::erf(a[i]);
  }
  return tmp;
}

inline simd::Float multiplysign(const simd::Float& x, const simd::Float& y) { // return x times sign of y
  return simd::Float(_mm_xor_ps(x._data, _mm_and_ps(hidden::SIGN_MASKf._data,y._data)));
}
  
inline simd::Float copysign(const simd::Float& x, const simd::Float& y) { // return abs(x) times sign of y
  return simd::Float(_mm_xor_ps(_mm_andnot_ps(hidden::SIGN_MASKf._data, x._data) , _mm_and_ps(hidden::SIGN_MASKf._data,y._data)));
}

inline simd::Float abs(const simd::Float& x) {
  return simd::Float(_mm_andnot_ps(hidden::SIGN_MASKf._data, x._data)); // !sign_mask & x
}

inline simd::Float max(const simd::Float& x, const simd::Float& y) {
  return simd::Float(_mm_max_ps(x._data,y._data));
}
  
inline simd::Float min(const simd::Float& x, const simd::Float& y) {
  return simd::Float(_mm_min_ps(x._data,y._data));
}

inline simd::Boolf isnan(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    if ( a[i] != a[i] ) {
      tmp[i] = 1.0;
    } else {
      tmp[i] = 0.0;
    }
  }
  return tmp > simd::Float(0.5f);
}

inline simd::Float if_then_else(const simd::Boolf& b, const simd::Float& v1, const simd::Float& v2) {
  return simd::Float( _mm_add_ps(_mm_and_ps(b._data,v1._data),_mm_andnot_ps(b._data,v2._data)) );
}

inline simd::Float if_then_else_zero(const simd::Boolf& b, const simd::Float& v) {
  return simd::Float( _mm_and_ps(b._data,v._data) );
}

} // namespace math
} // namespace stk
