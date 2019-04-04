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
static STK_MATH_FORCE_INLINE float CBRT(const float x) { return cbrtf(x); }
static const simd::Float SIGN_MASKf(-0.0);
static const simd::Float ZEROf(0.0);
}

STK_MATH_FORCE_INLINE simd::Float fmadd(const simd::Float& a, const simd::Float& b, const simd::Float& c) {
  return simd::Float(_mm512_add_ps(_mm512_mul_ps(a._data,b._data), c._data));
}

STK_MATH_FORCE_INLINE simd::Float sqrt(const simd::Float& x) {
  return simd::Float(_mm512_sqrt_ps(x._data));
}
  
STK_MATH_FORCE_INLINE simd::Float cbrt(const simd::Float& x) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm512_cbrt_ps(x._data));
#else
  simd::Float out;
  for (int i=0; i < simd::nfloats; ++i) {
    out[i] = hidden::CBRT(x[i]);
  }
  return out;
#endif
}

STK_MATH_FORCE_INLINE simd::Float log(const simd::Float& x) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm512_log_ps(x._data));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::log(x[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Float log10(const simd::Float& x) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::log10(x[n]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float exp(const simd::Float& x) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm512_exp_ps(x._data));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::exp(x[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, const simd::Float& y) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm512_pow_ps(x._data,y._data));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, float y) {
#if defined(__INTEL_COMPILER)
  return simd::Float(_mm512_pow_ps(x._data,_mm512_set1_ps(y)));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Float sin(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::sin(a[i]);
  }
  return tmp;
}
 
STK_MATH_FORCE_INLINE simd::Float cos(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::cos(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float tan(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::tan(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float sinh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::sinh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float cosh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::cosh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float tanh(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::tanh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float asin(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::asin(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float acos(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::acos(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float atan(const simd::Float& a) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::atan(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float atan2(const simd::Float& a, const simd::Float& b) {
  simd::Float tmp;
  for (int i=0; i < simd::nfloats; ++i) {
    tmp[i] = std::atan2(a[i], b[i]);
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

STK_MATH_FORCE_INLINE simd::Float multiplysign(const simd::Float& x, const simd::Float& y) { // return x times sign of y
  __m512i sign_y =  _mm512_and_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASKf._data), reinterpret_cast<const __m512i&>(y._data));
  __m512i sol = _mm512_xor_epi64(reinterpret_cast<const __m512i&>(x._data),sign_y);
  return simd::Float( reinterpret_cast<__m512&>(sol) );
}
  
STK_MATH_FORCE_INLINE simd::Float copysign(const simd::Float& x, const simd::Float& y) { // return abs(x) times sign of y
  __m512i pos_x = _mm512_andnot_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASKf._data), reinterpret_cast<const __m512i&>(x._data));
  __m512i sign_y = _mm512_and_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASKf._data), reinterpret_cast<const __m512i&>(y._data));
  __m512i sol = _mm512_xor_epi64(pos_x,sign_y);
  return simd::Float( reinterpret_cast<__m512&>(sol) );
}

STK_MATH_FORCE_INLINE simd::Float abs(const simd::Float& x) {  // !sign_mask & x
  __m512i sol = _mm512_andnot_epi64(  reinterpret_cast<const __m512i&>(hidden::SIGN_MASKf._data), reinterpret_cast<const __m512i&>(x._data) ); 
  return simd::Float( reinterpret_cast<__m512&>(sol) );
}

STK_MATH_FORCE_INLINE simd::Float min(const simd::Float& x, const simd::Float& y) {
  return simd::Float(_mm512_min_ps(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Float max(const simd::Float& x, const simd::Float& y) {
  return simd::Float(_mm512_max_ps(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Boolf isnan(const simd::Float& a) {
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

STK_MATH_FORCE_INLINE simd::Float if_then_else(const simd::Boolf& b, const simd::Float& v1, const simd::Float& v2) {
  return simd::Float( _mm512_mask_add_ps(v2._data, b._data,v1._data, hidden::ZEROf._data) );
}

STK_MATH_FORCE_INLINE simd::Float if_then_else_zero(const simd::Boolf& b, const simd::Float& v) {
  return simd::Float( _mm512_mask_add_ps(hidden::ZEROf._data, b._data, v._data, hidden::ZEROf._data) );
}

} // namespace math
} // namespace stk
