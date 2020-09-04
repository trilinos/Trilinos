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

#include "gtest/gtest.h"
#include <cmath>
#include "SimdFloatingPointFixture.hpp"

using SimdFloatMath = SimdFloatingPointFixture<stk::simd::Float, stk::simd::Float, stk::simd::nfloats>;

TEST_F(SimdFloatMath, copysign_posNeg)
{
  fill_constant(1.0, -2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](float a, float b){ return std::copysign(a, b); });
  verify();
}

TEST_F(SimdFloatMath, copysign_negPos)
{
  fill_constant(-2.0, 1.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](float a, float b){ return std::copysign(a, b); });
  verify();
}

TEST_F(SimdFloatMath, copysign_negNeg)
{
  fill_constant(-3.0, -4.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](float a, float b){ return std::copysign(a, b); });
  verify();
}

TEST_F(SimdFloatMath, copysign_posPos)
{
  fill_constant(6.0, 5.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](float a, float b){ return std::copysign(a, b); });
  verify();
}

TEST_F(SimdFloatMath, copysign_varying)
{
  if (is_scalar()) return;

  fill_varying_linear_plus_minus();
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](float a, float b){ return std::copysign(a, b); });
  verify();
}

TEST_F(SimdFloatMath, multiplysign_posNeg)
{
  fill_constant(1.0, -2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::multiplysign(a, b); });
  compute_expected_result([](float a, float b){ return a * std::copysign(1.0, b); });
  verify();
}

TEST_F(SimdFloatMath, multiplysign_negPos)
{
  fill_constant(-2.0, 1.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::multiplysign(a, b); });
  compute_expected_result([](float a, float b){ return a * std::copysign(1.0, b); });
  verify();
}

TEST_F(SimdFloatMath, multiplysign_negNeg)
{
  fill_constant(-3.0, -4.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::multiplysign(a, b); });
  compute_expected_result([](float a, float b){ return a * std::copysign(1.0, b); });
  verify();
}

TEST_F(SimdFloatMath, multiplysign_posPos)
{
  fill_constant(6.0, 5.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::multiplysign(a, b); });
  compute_expected_result([](float a, float b){ return a * std::copysign(1.0, b); });
  verify();
}

TEST_F(SimdFloatMath, multiplysign_varying)
{
  if (is_scalar()) return;

  fill_varying_linear_plus_minus();
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return stk::math::multiplysign(a, b); });
  compute_expected_result([](float a, float b){ return a * std::copysign(1.0, b); });
  verify();
}

using SimdFloatOperator = SimdFloatingPointFixture<stk::simd::Float, stk::simd::Boolf, stk::simd::nfloats>;

TEST_F(SimdFloatOperator, equal_aIsSmaller)
{
  fill_constant(1.0, 2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a == b; });
  compute_expected_result([](float a, float b){ return a == b; });
  verify();
}

TEST_F(SimdFloatOperator, equal_valuesAreEqual)
{
  fill_constant(2.0, 2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a == b; });
  compute_expected_result([](float a, float b){ return a == b; });
  verify();
}

TEST_F(SimdFloatOperator, equal_varying)
{
  if (is_scalar()) return;

  fill_varying_opposite_linear();
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a == b; });
  compute_expected_result([](float a, float b){ return a == b; });
  verify();
}

TEST_F(SimdFloatOperator, greaterThanEqual_aIsSmaller)
{
  fill_constant(1.0, 2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a >= b; });
  compute_expected_result([](float a, float b){ return a >= b; });
  verify();
}

TEST_F(SimdFloatOperator, greaterThanEqual_valuesAreEqual)
{
  fill_constant(2.0, 2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a >= b; });
  compute_expected_result([](float a, float b){ return a >= b; });
  verify();
}

TEST_F(SimdFloatOperator, greaterThanEqual_aIsLarger)
{
  fill_constant(3.0, -1.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a >= b; });
  compute_expected_result([](float a, float b){ return a >= b; });
  verify();
}

TEST_F(SimdFloatOperator, greaterThanEqual_varying)
{
  if (is_scalar()) return;

  fill_varying_opposite_linear();
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a >= b; });
  compute_expected_result([](float a, float b){ return a >= b; });
  verify();
}

TEST_F(SimdFloatOperator, lessThan_aIsSmaller)
{
  fill_constant(1.0, 2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a < b; });
  compute_expected_result([](float a, float b){ return a < b; });
  verify();
}

TEST_F(SimdFloatOperator, lessThan_valuesAreEqual)
{
  fill_constant(2.0, 2.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a < b; });
  compute_expected_result([](float a, float b){ return a < b; });
  verify();
}

TEST_F(SimdFloatOperator, lessThan_aIsLarger)
{
  fill_constant(3.0, -1.0);
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a < b; });
  compute_expected_result([](float a, float b){ return a < b; });
  verify();
}

TEST_F(SimdFloatOperator, lessThan_varying)
{
  if (is_scalar()) return;

  fill_varying_opposite_linear();
  compute_function([](const stk::simd::Float& a, const stk::simd::Float& b){ return a < b; });
  compute_expected_result([](float a, float b){ return a < b; });
  verify();
}

