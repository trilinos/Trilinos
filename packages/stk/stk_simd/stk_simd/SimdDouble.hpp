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

#ifndef STK_SIMD_DOUBLE_HPP
#define STK_SIMD_DOUBLE_HPP

namespace stk {
namespace simd {

struct Double {
  using native_simd_t = SIMD_NAMESPACE::native_simd<double>;
  
#if !defined(STK_VOLATILE_SIMD)
  static_assert(std::is_trivially_copyable<native_simd_t>::value,
                 "native_simd<double> must be trivially copyable");
  static_assert(std::is_trivially_destructible<native_simd_t>::value,
                 "native_simd<double> must be trivially destructible");
#endif

  STK_MATH_FORCE_INLINE Double() {}

  template <typename T, typename = typename std::enable_if<std::is_convertible<T, double>::value>::type>
  STK_MATH_FORCE_INLINE Double(const T x)
    : _data(static_cast<double>(x)) {
  }

  STK_MATH_FORCE_INLINE Double(const native_simd_t& x)
    : _data(x) {
  }

  STK_MATH_FORCE_INLINE Double(const Double& x)
    : _data(x._data) {
  }

#ifdef STK_VOLATILE_SIMD
  STK_MATH_FORCE_INLINE Double(const volatile Double& x)
    : _data(x._data) {
  }
#endif

  STK_MATH_FORCE_INLINE Double& operator= (const Double& x) {
    _data = x._data;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<std::is_convertible<T, double>::value>::type>
  STK_MATH_FORCE_INLINE Double& operator=(const T x) {
    _data = static_cast<double>(x);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const Double& a) {
    _data += a._data;
    return *this;
  }

#ifdef STK_VOLATILE_SIMD
  STK_MATH_FORCE_INLINE void operator+= (const volatile Double& a) volatile {
    _data.plus_equals(a._data);
  }
#endif

  STK_MATH_FORCE_INLINE Double& operator-= (const Double& a) {
    _data -= a._data;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const Double& a) {
    _data *= a._data;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const Double& a) {
    _data /= a._data;
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const double a) {
    _data += native_simd_t(a);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const double a) {
    _data -= native_simd_t(a);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const double a) {
    _data *= native_simd_t(a);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const double a) {
    _data /= native_simd_t(a);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double operator-() const {
    return Double(-_data);
  }

  STK_MATH_FORCE_INLINE double& operator[](int i) {
    assert(i >= 0 && static_cast<std::size_t>(i) < ndoubles);
    return (reinterpret_cast<double*>(&_data))[i];
  }

  STK_MATH_FORCE_INLINE const double& operator[](int i) const {
    assert(i >= 0 && static_cast<std::size_t>(i) < ndoubles);
    return (reinterpret_cast<const double*>(&_data))[i];
  }

  native_simd_t _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd
} // namespace stk

#endif //STK_SIMD_DOUBLE_HPP

