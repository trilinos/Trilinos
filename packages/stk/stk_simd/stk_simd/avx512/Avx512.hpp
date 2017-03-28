// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_AVX512_H
#define STK_SIMD_AVX512_H

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

#endif // STK_SIMD_AVX512_H
