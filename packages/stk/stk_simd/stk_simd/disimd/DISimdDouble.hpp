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

#ifndef STK_SIMD_DISIMD_DISIMDDOUBLE_HPP
#define STK_SIMD_DISIMD_DISIMDDOUBLE_HPP

namespace stk {
namespace simd {

struct Double {

  STK_MATH_FORCE_INLINE Double() {}


  template <typename T>
  STK_MATH_FORCE_INLINE Double(const T x, typename std::enable_if<std::is_convertible<T,double>::value, void*>::type=0)
    : _data(double(x)) {
  }

  STK_MATH_FORCE_INLINE Double(const SIMD_NAMESPACE::simd<double, SIMD_NAMESPACE::simd_abi::native>& x)
    : _data(x) {
  }

  STK_MATH_FORCE_INLINE Double(const Double& x)
    : _data(x._data) {
  }

  STK_MATH_FORCE_INLINE Double& operator= (const Double& x) {
    _data = x._data;
    return *this;
  }

  template <typename T>
  STK_MATH_FORCE_INLINE typename std::enable_if<std::is_convertible<T,double>::value, Double&>::type operator= (const T x) {
    _data = double(x);
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator+= (const Double& a) {
    _data += a._data;
    return *this;
  }

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
    _data = _data + a; //DI-QUESTION: operator+= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator-= (const double a) {
    _data = _data - a; //DI-QUESTION: operator-= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator*= (const double a) {
    _data = _data * a; //DI-QUESTION: operator*= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Double& operator/= (const double a) {
    _data = _data / a; //DI-QUESTION: operator/= ?
    return *this;
  }

  STK_MATH_FORCE_INLINE Double operator-() const {
    return (SIMD_NAMESPACE::simd<double, SIMD_NAMESPACE::simd_abi::native>(0.0) - _data);
  }

  STK_MATH_FORCE_INLINE double& operator[](int i) {return (reinterpret_cast<double*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const double& operator[](int i) const {return (reinterpret_cast<const double*>(&_data))[i];}
    
  STK_MATH_FORCE_INLINE int32_t& Int(int i) {return (reinterpret_cast<int32_t*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const int32_t& Int(int i) const {return (reinterpret_cast<const int32_t*>(&_data))[i];}

  STK_MATH_FORCE_INLINE uint32_t& UInt(int i) {return (reinterpret_cast<uint32_t*>(&_data))[i];}
  STK_MATH_FORCE_INLINE const uint32_t& UInt(int i) const {return (reinterpret_cast<const uint32_t*>(&_data))[i];}

  SIMD_NAMESPACE::simd<double, SIMD_NAMESPACE::simd_abi::native> _data; // the "_" means you should try not to use this directly
  // it is made public to avoid function call overhead 
  // and/or so the compiler doesn't have to use up one of
  // inlining depths (usually max inlining depth ~5)

};

} // namespace simd
} // namespace stk

#endif //STK_SIMD_DISIMD_DISIMDDOUBLE_HPP

