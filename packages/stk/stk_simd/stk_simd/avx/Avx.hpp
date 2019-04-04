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

inline double reduce_max(const Double& x) {
  double max = x[0];
  for (int i=1; i<4; ++i){
    max = max > x[i] ? max : x[i];
  }
  return max;
}

inline float reduce_max(const Float& x) {
  float max = x[0];
  for (int i=1; i<8; ++i){
    max = max > x[i] ? max : x[i];
  }
  return max;
}

inline double reduce_min(const Double& x) {
  double max = x[0];
  for (int i=1; i<4; ++i){
    max = max < x[i] ? max : x[i];
  }
  return max;
}

inline float reduce_min(const Float& x) {
  float max = x[0];
  for (int i=1; i<8; ++i){
    max = max < x[i] ? max : x[i];
  }
  return max;
}

}
}

#endif // STK_SIMD_AVX_H
