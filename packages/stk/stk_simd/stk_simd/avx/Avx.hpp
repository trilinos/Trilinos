// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_AVX_H
#define STK_SIMD_AVX_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <immintrin.h>
#include <stdio.h>
#include <cmath>
#include <assert.h>

namespace stk {
namespace simd {
constexpr int ndoubles = 4;
constexpr int nfloats = 8;
}
}

#include "./AvxDouble.hpp"
#include "./AvxFloat.hpp"
#include "./AvxBool.hpp"
#include "./AvxBoolF.hpp"

#include "./AvxDoubleOperators.hpp"
#include "./AvxDoubleLoadStore.hpp"
#include "./AvxDoubleMath.hpp"

#include "./AvxFloatOperators.hpp"
#include "./AvxFloatLoadStore.hpp"
#include "./AvxFloatMath.hpp"

namespace stk {
namespace simd {

inline double reduce_sum(const Double& x) {
  return x[0]+x[1]+x[2]+x[3];
}

inline float reduce_sum(const Float& x) {
  return x[0]+x[1]+x[2]+x[3]+
         x[4]+x[5]+x[6]+x[7];
}

}
}

#endif // STK_SIMD_AVX_H
