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

// IWYU pragma: private; include "Simd.hpp"

#ifndef STK_INCLUDE_ONLY_STK_SIMD_HEADER
static_assert(false, "Do not include simd impl files directly. Only include stk_simd/Simd.hpp");
#endif

#ifndef STK_SIMD_FLOAT_HPP
#define STK_SIMD_FLOAT_HPP

namespace stk {
namespace simd {

struct Float {
  using native_simd_t = SIMD_NAMESPACE::native_simd<float>;

  STK_MATH_FORCE_INLINE Float() {}

  template <typename T>
  STK_MATH_FORCE_INLINE Float(const T x, typename std::enable_if<std::is_convertible<T,float>::value, void*>::type=0)
    : _data(static_cast<float>(x)) {
  }

  STK_MATH_FORCE_INLINE Float(const native_simd_t& x)
    : _data(x.get()) {
  }

  STK_MATH_FORCE_INLINE Float(const Float& x)
    : _data(x._data.get()) {
  }

  STK_MATH_FORCE_INLINE Float& operator= (const Float& x) {
    _data = x._data;
    return *this;
  }

  template <typename T>
  STK_MATH_FORCE_INLINE typename std::enable_if<std::is_convertible<T,float>::value, Float&>::type operator= (const T x) {
    _data = static_cast<float>(x);
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator+= (const Float& a) {
    _data += a._data;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator-= (const Float& a) {
    _data -= a._data;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator*= (const Float& a) {
    _data *= a._data;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator/= (const Float& a) {
    _data /= a._data;
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator+= (const float a) {
    _data = _data + a; //DI-QUESTION: operator+= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator-= (const float a) {
    _data = _data - a; //DI-QUESTION: operator-= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator*= (const float a) {
    _data = _data * a; //DI-QUESTION: operator*= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Float& operator/= (const float a) {
    _data = _data / a; //DI-QUESTION: operator/= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Float operator-() const {
    return - _data;
  }

  STK_MATH_FORCE_INLINE float& operator[](int i) {return (reinterpret_cast<float*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const float& operator[](int i) const {return (reinterpret_cast<const float*>(&_data))[i];}

  native_simd_t _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd
} // namespace stk

#endif //STK_SIMD_FLOAT_HPP

