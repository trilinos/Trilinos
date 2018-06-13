// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_FW_H
#define STK_SIMD_FW_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <stdio.h>
#include <cmath>
#include <assert.h>
#include "NoSimdSizes.hpp"

#include "NoSimdDouble.hpp"
#include "NoSimdFloat.hpp"
#include "NoSimdBool.hpp"
#include "NoSimdBoolF.hpp"

#include "NoSimdDoubleOperators.hpp"
#include "NoSimdDoubleLoadStore.hpp"
#include "NoSimdDoubleMath.hpp"

#include "NoSimdFloatOperators.hpp"
#include "NoSimdFloatLoadStore.hpp"
#include "NoSimdFloatMath.hpp"

namespace stk {
namespace simd {

STK_MATH_FORCE_INLINE double reduce_sum(const Double& x) {
  return x[0];
}

STK_MATH_FORCE_INLINE float reduce_sum(const Float& x) {
  return x[0]+x[1];
}

}
}

#endif // STK_SIMD_FW_H
