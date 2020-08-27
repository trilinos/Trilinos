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
// 

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <array>
#include "./NoSimdSizes.hpp"

namespace stk {
namespace simd {

struct Double {

  STK_MATH_FORCE_INLINE Double() {}

  template <typename T>
  STK_MATH_FORCE_INLINE Double(const T x, typename std::enable_if<std::is_convertible<T,double>::value, void*>::type=0) {

    for (int i=0; i < ndoubles; ++i) _data[i] = double(x);
  }

  STK_MATH_FORCE_INLINE Double(const Double& x) {

    for (int i=0; i < ndoubles; ++i) _data[i] = x._data[i];  
  }

  STK_MATH_FORCE_INLINE Double& operator= (const Double& x) {

    for (int i=0; i < ndoubles; ++i) _data[i] = x._data[i];
    return *this;
  }

  template <typename T>
  STK_MATH_FORCE_INLINE typename std::enable_if<std::is_convertible<T,double>::value, Double&>::type operator= (const T x) {

    for (int i=0; i < ndoubles; ++i) _data[i] = double(x);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] += a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] -= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] *= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const Double& a) {

    for (int i=0; i < ndoubles; ++i) _data[i] /= a[i];
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] += a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] -= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] *= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const double a) {

    for (int i=0; i < ndoubles; ++i) _data[i] /= a;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double operator-() const {

    Double tmp;
    for (int i=0; i < ndoubles; ++i) tmp[i] = -_data[i];
    return tmp;
  }

  STK_MATH_FORCE_INLINE double& operator[](int i) {return _data[i];}
  STK_MATH_FORCE_INLINE const double& operator[](int i) const {return _data[i];}

  double _data[simd::ndoubles];
};

} // namespace simd
} // namespace stk

