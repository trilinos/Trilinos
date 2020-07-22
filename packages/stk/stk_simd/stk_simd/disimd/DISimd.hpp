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

#ifndef STK_SIMD_DISIMD_HPP
#define STK_SIMD_DISIMD_HPP

#ifndef SIMD_ALWAYS_INLINE
//currently necessary to avoid the 'always_inline' defined in simd.hpp
#define SIMD_ALWAYS_INLINE __attribute__((always_inline))
#endif

#ifdef STK_SIMD_NONE
#define native native_junk
#define native_simd native_simd_junk
#endif

#include "./simd.hpp"

#ifdef STK_SIMD_NONE
#undef native
#undef native_simd

namespace SIMD_NAMESPACE {

namespace simd_abi {
using native = scalar;
}

template <class T>
using native_simd = simd<T, simd_abi::native>;

}
#endif

namespace stk {
namespace simd {
constexpr int ndoubles = SIMD_NAMESPACE::simd<double, SIMD_NAMESPACE::simd_abi::native>::size();
constexpr int nfloats = SIMD_NAMESPACE::simd<float, SIMD_NAMESPACE::simd_abi::native>::size();
}
}

#include "./DISimdDouble.hpp"
#include "./DISimdFloat.hpp"
#include "./DISimdBool.hpp"
#include "./DISimdBoolF.hpp"
//
#include "./DISimdDoubleOperators.hpp"
#include "./DISimdDoubleLoadStore.hpp"
#include "./DISimdDoubleMath.hpp"
//
#include "./DISimdFloatOperators.hpp"
#include "./DISimdFloatLoadStore.hpp"
#include "./DISimdFloatMath.hpp"

namespace stk {
namespace simd {

STK_MATH_FORCE_INLINE double reduce_sum(const Double& x) {
  double sum = x[0];
  for (int i=1; i<ndoubles; ++i) {
    sum += x[i];
  }
  return sum;
}

STK_MATH_FORCE_INLINE float reduce_sum(const Float& x) {
  double sum = x[0];
  for (int i=1; i<nfloats; ++i) {
    sum += x[i];
  }
  return sum;
}

STK_MATH_FORCE_INLINE double reduce_max(const Double& x) {
  double max = x[0];
  for (int i=1; i<ndoubles; ++i){
    max = max > x[i] ? max : x[i];
  }
  return max;
}

STK_MATH_FORCE_INLINE float reduce_max(const Float& x) {
  float max = x[0];
  for (int i=1; i<nfloats; ++i){
    max = max > x[i] ? max : x[i];
  }
  return max;
}

STK_MATH_FORCE_INLINE double reduce_min(const Double& x) {
  double max = x[0];
  for (int i=1; i<ndoubles; ++i){
    max = max < x[i] ? max : x[i];
  }
  return max;
}

STK_MATH_FORCE_INLINE float reduce_min(const Float& x) {
  float max = x[0];
  for (int i=1; i<nfloats; ++i){
    max = max < x[i] ? max : x[i];
  }
  return max;
}

}
}

#endif // STK_SIMD_DISIMD_HPP
