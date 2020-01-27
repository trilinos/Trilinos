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

#ifndef STK_SIMD_SSE_H
#define STK_SIMD_SSE_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <emmintrin.h>
#include <stdio.h>
#include <cmath>
#include <assert.h>

namespace stk {
namespace simd {
constexpr int ndoubles = 2;
constexpr int nfloats  = 4;
}
}

#include "./SseDouble.hpp"
#include "./SseFloat.hpp"
#include "./SseBool.hpp"
#include "./SseBoolF.hpp"

#include "./SseDoubleOperators.hpp"
#include "./SseDoubleLoadStore.hpp"
#include "./SseDoubleMath.hpp"

#include "./SseFloatOperators.hpp"
#include "./SseFloatLoadStore.hpp"
#include "./SseFloatMath.hpp"

namespace stk {
namespace simd {

inline double reduce_sum(const Double& x) {
  return x[0]+x[1];
}

inline float reduce_sum(const Float& x) {
  return x[0]+x[1]+x[2]+x[3];
}

inline double reduce_max(const Double& x) {
  double max = x[0];
  for (int i=1; i<2; ++i){
    max = max > x[i] ? max : x[i];
  }
  return max;
}

inline float reduce_max(const Float& x) {
  float max = x[0];
  for (int i=1; i<4; ++i){
    max = max > x[i] ? max : x[i];
  }
  return max;
}

inline double reduce_min(const Double& x) {
  double max = x[0];
  for (int i=1; i<2; ++i){
    max = max < x[i] ? max : x[i];
  }
  return max;
}

inline float reduce_min(const Float& x) {
  float max = x[0];
  for (int i=1; i<4; ++i){
    max = max < x[i] ? max : x[i];
  }
  return max;
}

}
}

#endif // STK_SIMD_SSE_H
