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

#ifndef STK_SIMD_AVX512_H
#define STK_SIMD_AVX512_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <immintrin.h>
#include <stdio.h>
#include <cmath>
#include <assert.h>

namespace stk {
namespace simd {
constexpr int ndoubles = 8;
constexpr int nfloats = 16;
}
}

#include "./Avx512Double.hpp"
#include "./Avx512Float.hpp"
#include "./Avx512Bool.hpp"
#include "./Avx512BoolF.hpp"

#include "./Avx512DoubleOperators.hpp"
#include "./Avx512DoubleLoadStore.hpp"
#include "./Avx512DoubleMath.hpp"

#include "./Avx512FloatOperators.hpp"
#include "./Avx512FloatLoadStore.hpp"
#include "./Avx512FloatMath.hpp"

namespace stk {
namespace simd {

inline double reduce_sum(const Double& x) {
  return _mm512_reduce_add_pd(x._data);
}

inline float reduce_sum(const Float& x) {
  return _mm512_reduce_add_ps(x._data);
}

inline double reduce_max(const Double& x) {
  return _mm512_reduce_max_pd(x._data);
}

inline float reduce_max(const Float& x) {
  return _mm512_reduce_max_ps(x._data);
}

inline double reduce_min(const Double& x) {
  return _mm512_reduce_min_pd(x._data);
}

inline float reduce_min(const Float& x) {
  return _mm512_reduce_min_ps(x._data);
}

}
}

#endif // STK_SIMD_AVX512_H
