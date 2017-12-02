// Copyright 2013 Sandia Corporation, Albuquerque, NM.

#ifndef STK_SIMD_FW_H
#define STK_SIMD_FW_H

// IWYU pragma: private, include <stk_simd/Simd.hpp>

#include <stdio.h>
#include <cmath>
#include <assert.h>
#include "./FwSizes.hpp"

#include "./FwDouble.hpp"
#include "./FwFloat.hpp"
#include "./FwBool.hpp"
#include "./FwBoolF.hpp"

#include "./FwDoubleOperators.hpp"
#include "./FwDoubleLoadStore.hpp"
#include "./FwDoubleMath.hpp"

#include "./FwFloatOperators.hpp"
#include "./FwFloatLoadStore.hpp"
#include "./FwFloatMath.hpp"

namespace stk {
namespace simd {

STK_MATH_FORCE_INLINE double reduce_sum(const Double& x) {
  return x[0]+x[1]+x[2]+x[3];
}

STK_MATH_FORCE_INLINE float reduce_sum(const Float& x) {
  return x[0]+x[1]+x[2]+x[3]+
         x[4]+x[5]+x[6]+x[7];
}

}
}

#endif // STK_SIMD_FW_H
