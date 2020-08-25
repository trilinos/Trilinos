// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#ifndef STK_SIMD_DISIMD_DISIMDDOUBLEMATH_HPP
#define STK_SIMD_DISIMD_DISIMDDOUBLEMATH_HPP

namespace stk {
namespace math {

namespace hidden {
STK_MATH_FORCE_INLINE static double CBRT(const double x) { return cbrt(x); }
static const simd::Double SIGN_MASK(-0.0);
}

STK_MATH_FORCE_INLINE simd::Double fmadd(const simd::Double& a, const simd::Double& b, const simd::Double& c) {
  return simd::Double(fma(a._data, b._data, c._data));
}

//STK_MATH_FORCE_INLINE simd::Double fmadd(const double a, const simd::Double& b, const simd::Double& c) {
//#if defined(__AVX2__)
//  return simd::Double(_mm256_fmadd_pd(_mm256_set1_pd(a), b._data, c._data));
//#else
//  return simd::Double(
//    _mm256_add_pd(
//      _mm256_mul_pd(_mm256_set1_pd(a), b._data), c._data)
//  );
//#endif
//}
//
//STK_MATH_FORCE_INLINE simd::Double fmadd(const double a, const double b, const simd::Double& c) {
//#if defined(__AVX2__)
//  return simd::Double(_mm256_fmadd_pd(_mm256_set1_pd(a), _mm256_set1_pd(b), c._data));
//#else 
//  return simd::Double(
//    _mm256_add_pd(
//      _mm256_mul_pd(_mm256_set1_pd(a), _mm256_set1_pd(b)), c._data)
//  );
//#endif
//}
//
//STK_MATH_FORCE_INLINE simd::Double fmadd(const double a, const simd::Double& b, const double c) {
//#if defined(__AVX2__)
//  return simd::Double(_mm256_fmadd_pd(_mm256_set1_pd(a), b._data, _mm256_set1_pd(c)));
//#else
//  return simd::Double(
//    _mm256_add_pd(
//      _mm256_mul_pd(_mm256_set1_pd(a), b._data), _mm256_set1_pd(c))
//  );
//#endif
//}
//
//STK_MATH_FORCE_INLINE simd::Double fmadd(const simd::Double& a, const double b, const simd::Double& c) {
//#if defined(__AVX2__)
//  return simd::Double(_mm256_fmadd_pd(a._data, _mm256_set1_pd(b), c._data));
//#else
//  return simd::Double(
//    _mm256_add_pd(
//      _mm256_mul_pd(a._data, _mm256_set1_pd(b)), c._data)
//  );
//#endif
//}
//
//STK_MATH_FORCE_INLINE simd::Double fmadd(const simd::Double& a, const double b, const double c) {
//#if defined(__AVX2__)
//  return simd::Double(_mm256_fmadd_pd(a._data, _mm256_set1_pd(b), _mm256_set1_pd(c)));
//#else
//  return simd::Double(
//    _mm256_add_pd(
//      _mm256_mul_pd(a._data, _mm256_set1_pd(b)), _mm256_set1_pd(c))
//  );
//#endif
//}
//
//STK_MATH_FORCE_INLINE simd::Double fmadd(const simd::Double& a, const simd::Double& b, const double c) {
//#if defined(__AVX2__)
//  return simd::Double(_mm256_fmadd_pd(a._data, b._data, _mm256_set1_pd(c)));
//#else
//  return simd::Double(
//    _mm256_add_pd(
//      _mm256_mul_pd(a._data, b._data), _mm256_set1_pd(c))
//  );
//#endif
//}

STK_MATH_FORCE_INLINE simd::Double sqrt(const simd::Double& x) {
  return simd::Double(SIMD_NAMESPACE::sqrt(x._data));
}
  
STK_MATH_FORCE_INLINE simd::Double cbrt(const simd::Double& x) {
#if defined(__INTEL_COMPILER) && !defined(STK_SIMD_NONE)
  return simd::Double(SIMD_NAMESPACE::cbrt(x._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = hidden::CBRT(x[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Double log(const simd::Double& x) {
#if defined(__INTEL_COMPILER) && !defined(STK_SIMD_NONE)
  return simd::Double(SIMD_NAMESPACE::log(x._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::log(x[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Double log10(const simd::Double& x) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::log10(x[n]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double exp(const simd::Double& x) {
#if defined(__INTEL_COMPILER) && !defined(STK_SIMD_NONE)
  return simd::Double(SIMD_NAMESPACE::exp(x._data));
#else
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::exp(x[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const simd::Double& y) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n],y[n]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const double y) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double pow(const simd::Double& x, const int y) {
  simd::Double tmp;
  for (int n=0; n < simd::ndoubles; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double sin(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::sin(a[i]);
  }
  return tmp;
}
 
STK_MATH_FORCE_INLINE simd::Double cos(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::cos(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double tan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::tan(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double sinh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::sinh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double cosh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::cosh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double tanh(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::tanh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double asin(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::asin(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double acos(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::acos(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double atan(const simd::Double& a) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::atan(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Double atan2(const simd::Double& a, const simd::Double& b) {
  simd::Double tmp;
  for (int i=0; i < simd::ndoubles; ++i) {
    tmp[i] = std::atan2(a[i],b[i]);
  }
  return tmp;
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
  return simd::Double(SIMD_NAMESPACE::multiplysign(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double copysign(const simd::Double& x, const simd::Double& y) { // return abs(x) times sign of y
  return simd::Double(SIMD_NAMESPACE::copysign(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double abs(const simd::Double& x) {
  return simd::Double(SIMD_NAMESPACE::abs(x._data));
}

STK_MATH_FORCE_INLINE simd::Double min(const simd::Double& x, const simd::Double& y) {
  return simd::Double(SIMD_NAMESPACE::min(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Double max(const simd::Double& x, const simd::Double& y) {
  return simd::Double(SIMD_NAMESPACE::max(x._data, y._data));
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
  return tmp > simd::Double(0.5);
}

STK_MATH_FORCE_INLINE simd::Double if_then_else(const simd::Bool& b, const simd::Double& v1, const simd::Double& v2) {
  return simd::Double(SIMD_NAMESPACE::choose(b._data, v1._data, v2._data));
}

STK_MATH_FORCE_INLINE simd::Double if_then_else_zero(const simd::Bool& b, const simd::Double& v) {
  return simd::Double(SIMD_NAMESPACE::choose(b._data, v._data, SIMD_NAMESPACE::simd<double, SIMD_NAMESPACE::simd_abi::native>(0.0)));
}

} // namespace math
} // namespace stk

#endif //STK_SIMD_DISIMD_DISIMDDOUBLEMATH_HPP

