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

// IWYU pragma: private; include "Simd.hpp"

#ifndef STK_INCLUDE_ONLY_STK_SIMD_HEADER
static_assert(false, "Do not include simd impl files directly. Only include stk_simd/Simd.hpp");
#endif

#ifndef STK_SIMD_FLOATMATH_HPP
#define STK_SIMD_FLOATMATH_HPP

namespace stk {
namespace math {

namespace hidden {
STK_MATH_FORCE_INLINE static float CBRT(const float x) { return cbrt(x); }
static const simd::Float SIGN_MASKF(-0.0);
}

STK_MATH_FORCE_INLINE simd::Float fmadd(const simd::Float& a, const simd::Float& b, const simd::Float& c) {
  return simd::Float(fma(a._data, b._data, c._data));
}

STK_MATH_FORCE_INLINE simd::Float sqrt(const simd::Float& x) {
  return simd::Float(SIMD_NAMESPACE::sqrt(x._data));
}
  
STK_MATH_FORCE_INLINE simd::Float cbrt(const simd::Float& x) {
#if (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) && !defined(USE_STK_SIMD_NONE)
  return simd::Float(SIMD_NAMESPACE::cbrt(x._data));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = hidden::CBRT(x[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Float log(const simd::Float& x) {
#if (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) && !defined(USE_STK_SIMD_NONE)
  return simd::Float(SIMD_NAMESPACE::log(x._data));
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
#if (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)) && !defined(USE_STK_SIMD_NONE)
  return simd::Float(SIMD_NAMESPACE::exp(x._data));
#else
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::exp(x[n]);
  }
  return tmp;
#endif
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, const simd::Float& y) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y[n]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, const float y) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float pow(const simd::Float& x, const int y) {
  simd::Float tmp;
  for (int n=0; n < simd::nfloats; ++n) {
    tmp[n] = std::pow(x[n],y);
  }
  return tmp;
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
    tmp[i] = std::atan2(a[i],b[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float asinh(const simd::Float &a) {
  simd::Float tmp;
  for (int i = 0; i < simd::nfloats; ++i) {
    tmp[i] = std::asinh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float acosh(const simd::Float &a) {
  simd::Float tmp;
  for (int i = 0; i < simd::nfloats; ++i) {
    tmp[i] = std::acosh(a[i]);
  }
  return tmp;
}

STK_MATH_FORCE_INLINE simd::Float atanh(const simd::Float &a) {
  simd::Float tmp;
  for (int i = 0; i < simd::nfloats; ++i) {
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
  return simd::Float(SIMD_NAMESPACE::multiplysign(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Float copysign(const simd::Float& x, const simd::Float& y) { // return abs(x) times sign of y
  return simd::Float(SIMD_NAMESPACE::copysign(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Float abs(const simd::Float& x) {
  return simd::Float(SIMD_NAMESPACE::abs(x._data));
}

STK_MATH_FORCE_INLINE simd::Float min(const simd::Float& x, const simd::Float& y) {
  return simd::Float(SIMD_NAMESPACE::min(x._data, y._data));
}

STK_MATH_FORCE_INLINE simd::Float max(const simd::Float& x, const simd::Float& y) {
  return simd::Float(SIMD_NAMESPACE::max(x._data, y._data));
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
  return tmp > simd::Float(0.5);
}

STK_MATH_FORCE_INLINE simd::Float if_then_else(const simd::Boolf& b, const simd::Float& v1, const simd::Float& v2) {
  return simd::Float(SIMD_NAMESPACE::choose(b._data, v1._data, v2._data));
}

STK_MATH_FORCE_INLINE simd::Float if_then_else_zero(const simd::Boolf& b, const simd::Float& v) {
  return simd::Float(SIMD_NAMESPACE::choose(b._data, v._data, SIMD_NAMESPACE::native_simd<float>(0.0)));
}

} // namespace math
} // namespace stk

#endif //STK_SIMD_FLOATMATH_HPP
