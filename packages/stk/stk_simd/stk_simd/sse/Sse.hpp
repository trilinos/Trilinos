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

#endif // STK_SIMD_SSE_H
