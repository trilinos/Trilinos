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

#include <cmath>
#include <gtest/gtest.h>

#include "SimdFloatingPointFixture.hpp"
#include "SimdNameGenerator.hpp"

template <typename T>
using SimdFloatingPointMath = SimdFloatingPointFixture<T, T>;
using FloatingPointTypes = ::testing::Types<stk::simd::Double, stk::simd::Float>;
TYPED_TEST_SUITE(SimdFloatingPointMath, FloatingPointTypes, SimdNameGenerator);

TYPED_TEST(SimdFloatingPointMath, copysign_posNeg)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(1.0, -2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::copysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return std::copysign(a, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, copysign_negPos)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(-2.0, 1.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::copysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return std::copysign(a, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, copysign_negNeg)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(-3.0, -4.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::copysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return std::copysign(a, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, copysign_posPos)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(6.0, 5.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::copysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return std::copysign(a, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, copysign_varying)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  if (this->is_scalar()) return;

  this->fill_varying_linear_plus_minus();
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::copysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return std::copysign(a, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, multiplysign_posNeg)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(1.0, -2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::multiplysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a * std::copysign(1.0, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, multiplysign_negPos)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(-2.0, 1.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::multiplysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a * std::copysign(1.0, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, multiplysign_negNeg)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(-3.0, -4.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::multiplysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a * std::copysign(1.0, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, multiplysign_posPos)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(6.0, 5.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::multiplysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a * std::copysign(1.0, b); });
  this->verify();
}

TYPED_TEST(SimdFloatingPointMath, multiplysign_varying)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  if (this->is_scalar()) return;

  this->fill_varying_linear_plus_minus();
  this->compute_function([](const SimdType& a, const SimdType& b){ return stk::math::multiplysign(a, b); });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a * std::copysign(1.0, b); });
  this->verify();
}


template <typename T>
using SimdFloatingPointOperator = SimdFloatingPointFixture<T, typename stk::BoolT<T>::type>;
using FloatingPointTypes = ::testing::Types<stk::simd::Double, stk::simd::Float>;
TYPED_TEST_SUITE(SimdFloatingPointOperator, FloatingPointTypes,);

TYPED_TEST(SimdFloatingPointOperator, equal_aIsSmaller)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(1.0, 2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a == b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a == b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, equal_valuesAreEqual)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(2.0, 2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a == b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a == b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, equal_varying)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  if (this->is_scalar()) return;

  this->fill_varying_opposite_linear();
  this->compute_function([](const SimdType& a, const SimdType& b){ return a == b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a == b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, greaterThanEqual_aIsSmaller)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(1.0, 2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a >= b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a >= b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, greaterThanEqual_valuesAreEqual)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(2.0, 2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a >= b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a >= b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, greaterThanEqual_aIsLarger)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(3.0, -1.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a >= b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a >= b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, greaterThanEqual_varying)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  if (this->is_scalar()) return;

  this->fill_varying_opposite_linear();
  this->compute_function([](const SimdType& a, const SimdType& b){ return a >= b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a >= b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, lessThan_aIsSmaller)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(1.0, 2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a < b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a < b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, lessThan_valuesAreEqual)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(2.0, 2.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a < b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a < b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, lessThan_aIsLarger)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  this->fill_constant(3.0, -1.0);
  this->compute_function([](const SimdType& a, const SimdType& b){ return a < b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a < b; });
  this->verify();
}

TYPED_TEST(SimdFloatingPointOperator, lessThan_varying)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;

  if (this->is_scalar()) return;

  this->fill_varying_opposite_linear();
  this->compute_function([](const SimdType& a, const SimdType& b){ return a < b; });
  this->compute_expected_result([](ScalarType a, ScalarType b){ return a < b; });
  this->verify();
}
