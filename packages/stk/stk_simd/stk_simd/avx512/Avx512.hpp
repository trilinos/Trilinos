// Copyright 2013 Sandia Corporation, Albuquerque, NM.

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

}
}

#endif // STK_SIMD_AVX512_H
