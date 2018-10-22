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
STK_MATH_FORCE_INLINE static double CBRT(const double x) { return cbrt(x); }
static const simd::Double SIGN_MASK(-0.0);
static const simd::Double ZERO(0.0);
}

STK_MATH_FORCE_INLINE simd::Double fmadd(const simd::Double& a, const simd::Double& b, const simd::Double& c) {
  return simd::Double(_mm512_fmadd_pd(a._data, b._data, c._data));
}

STK_MATH_FORCE_INLINE simd::Double sqrt(const simd::Double& x) {
  return simd::Double(_mm512_sqrt_pd(x._data));
}
  
STK_MATH_FORCE_INLINE simd::Double cbrt(const simd::Double& x) {
  return simd::Double(_mm512_cbrt_pd(x._data));
}

STK_MATH_FORCE_INLINE simd::Double log(const simd::Double& x) {
  return simd::Double(_mm512_log_pd(x._data));
}

STK_MATH_FORCE_INLINE simd::Double log10(const simd::Double& x) {
  return simd::Double(_mm512_log10_pd(x._data));
}

STK_MATH_FORCE_INLINE simd::Double exp(const simd::Double& x) {
  return simd::Double(_mm512_exp_pd(x._data));
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const simd::Double& y) {
  return simd::Double(_mm512_pow_pd(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, double y) {
  return simd::Double(_mm512_pow_pd(x._data, _mm512_set1_pd(y)));
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, int y) {
  if (y==1) return x;
  if (y==2) return x*x;
  if (y==3) return x*x*x;
  return pow(x, double(y));
}

STK_MATH_FORCE_INLINE simd::Double sin(const simd::Double& a) {
  return simd::Double(_mm512_sin_pd(a._data));
}
 
STK_MATH_FORCE_INLINE simd::Double cos(const simd::Double& a) {
  return simd::Double(_mm512_cos_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double tan(const simd::Double& a) {
  return simd::Double(_mm512_tan_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double sinh(const simd::Double& a) {
  return simd::Double(_mm512_sinh_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double cosh(const simd::Double& a) {
  return simd::Double(_mm512_cosh_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double tanh(const simd::Double& a) {
  return simd::Double(_mm512_tanh_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double asin(const simd::Double& a) {
  return simd::Double(_mm512_asin_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double acos(const simd::Double& a) {
  return simd::Double(_mm512_acos_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double atan(const simd::Double& a) {
  return simd::Double(_mm512_atan_pd(a._data));
}

STK_MATH_FORCE_INLINE simd::Double atan2(const simd::Double& a, const simd::Double& b) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::atan2(a[i],b[i]);
  }
  return tmp;
  // the intrinsic is slower last time it was tested on 4/28/17
  //return simd::Double(_mm512_atan2_pd(a._data, b._data));
}

STK_MATH_FORCE_INLINE simd::Double asinh(const simd::Double &a) {
  simd::Double tmp;
  for (int i = 0; i < simd::ndoubles; ++i) {
    tmp[i] = std::asinh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double acosh(const simd::Double &a) {
  simd::Double tmp;
  for (int i = 0; i < simd::ndoubles; ++i) {
    tmp[i] = std::acosh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double atanh(const simd::Double &a) {
  simd::Double tmp;
  for (int i = 0; i < simd::ndoubles; ++i) {
    tmp[i] = std::atanh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double erf(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::erf(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double multiplysign(const simd::Double& x, const simd::Double& y) { // return x times sign of y
  __m512i sign_y =  _mm512_and_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASK._data),reinterpret_cast<const __m512i&>(y._data));
  __m512i sol = _mm512_xor_epi64(reinterpret_cast<const __m512i&>(x._data),sign_y);
  return simd::Double( reinterpret_cast<__m512d&>(sol) );
}
  
STK_MATH_FORCE_INLINE simd::Double copysign(const simd::Double& x, const simd::Double& y) { // return abs(x) times sign of y
  __m512i pos_x = _mm512_andnot_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASK._data), reinterpret_cast<const __m512i&>(x._data));
  __m512i sign_y = _mm512_and_epi64(reinterpret_cast<const __m512i&>(hidden::SIGN_MASK._data),reinterpret_cast<const __m512i&>(y._data));
  __m512i sol = _mm512_xor_epi64(pos_x,sign_y);
  return simd::Double( reinterpret_cast<__m512d&>(sol) );
}

STK_MATH_FORCE_INLINE simd::Double abs(const simd::Double& x) {
  return simd::Double(_mm512_abs_pd(x._data));
}

STK_MATH_FORCE_INLINE simd::Double min(const simd::Double& x, const simd::Double& y) {
  return simd::Double(_mm512_min_pd(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double max(const simd::Double& x, const simd::Double& y) {
  return simd::Double(_mm512_max_pd(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Bool isnan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    if ( a[i] != a[i] ) {
      tmp[i] = 1.0;
    } else {
      tmp[i] = 0.0;
    }
  }
  return tmp > simd::Double(0.5f);
}

STK_MATH_FORCE_INLINE simd::Double if_then_else(const simd::Bool& b, const simd::Double& v1, const simd::Double& v2) {
  return simd::Double( _mm512_mask_add_pd(v2._data, b._data, v1._data, hidden::ZERO._data) );
}

STK_MATH_FORCE_INLINE simd::Double if_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  return simd::Double( _mm512_mask_add_pd(hidden::ZERO._data, b._data, v._data, hidden::ZERO._data) );
}

} // namespace math
} // namespace stk
