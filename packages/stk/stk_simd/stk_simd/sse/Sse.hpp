// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_SSE_H
#define STK_SIMD_SSE_H

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

// IWYU pragma: begin_exports
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
// IWYU pragma: end_exports

namespace stk {
namespace simd {

inline double reduce_sum(const Double& x) {
  return x[0]+x[1];
}

inline float reduce_sum(const Float& x) {
  return x[0]+x[1]+x[2]+x[3];
}

}
}

#endif // STK_SIMD_SSE_H
