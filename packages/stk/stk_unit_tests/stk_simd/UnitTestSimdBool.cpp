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
#include "SimdBoolFixture.hpp"


using BoolTypes = ::testing::Types<stk::simd::Bool, stk::simd::Boolf>;
TYPED_TEST_SUITE(SimdBoolFixture, BoolTypes,);

TYPED_TEST(SimdBoolFixture, selectByLane_allTrue)
{
  this->fill_constant(true);
  this->use_simd_value_as_result();
  this->set_expected_result([](int /*i*/){ return true; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, selectByLane_allFalse)
{
  this->fill_constant(false);
  this->use_simd_value_as_result();
  this->set_expected_result([](int /*i*/){ return false; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, selectByLane_someTrue)
{
  if (this->is_scalar()) return;

  this->fill_alternating();
  this->use_simd_value_as_result();
  this->set_expected_result([](int i){ return (i%2 == 0); });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorOr_TrueTrue)
{
  this->fill_constant(true, true);
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a || b; });
  this->compute_expected_result([](bool a, bool b){ return a || b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorOr_FalseFalse)
{
  this->fill_constant(false, false);
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a || b; });
  this->compute_expected_result([](bool a, bool b){ return a || b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorOr_TrueFalse)
{
  this->fill_constant(true, false);
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a || b; });
  this->compute_expected_result([](bool a, bool b){ return a || b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorOr_Alternating)
{
  if (this->is_scalar()) return;

  this->fill_alternating();
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a || b; });
  this->compute_expected_result([](bool a, bool b){ return a || b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorAnd_TrueTrue)
{
  this->fill_constant(true, true);
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a && b; });
  this->compute_expected_result([](bool a, bool b){ return a && b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorAnd_FalseFalse)
{
  this->fill_constant(false, false);
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a && b; });
  this->compute_expected_result([](bool a, bool b){ return a && b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorAnd_TrueFalse)
{
  this->fill_constant(true, false);
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a && b; });
  this->compute_expected_result([](bool a, bool b){ return a && b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorAnd_Alternating)
{
  if (this->is_scalar()) return;

  this->fill_alternating();
  this->compute_function([](const TypeParam& a, const TypeParam& b){ return a && b; });
  this->compute_expected_result([](bool a, bool b){ return a && b; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorNot_allTrue)
{
  this->fill_constant(true);
  this->compute_function([](const TypeParam& a){ return !a; });
  this->compute_expected_result([](bool a){ return !a; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorNot_allFalse)
{
  this->fill_constant(false);
  this->compute_function([](const TypeParam& a){ return !a; });
  this->compute_expected_result([](bool a){ return !a; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, operatorNot_someTrue)
{
  if (this->is_scalar()) return;

  this->fill_alternating();
  this->compute_function([](const TypeParam& a){ return !a; });
  this->compute_expected_result([](bool a){ return !a; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, allOf_allTrue)
{
  this->fill_constant(true);
  this->compute_function([](const TypeParam& a)->bool{ return stk::simd::are_all(a); });
  this->set_expected_result([](){ return true; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, allOf_allFalse)
{
  this->fill_constant(false);
  this->compute_function([](const TypeParam& a)->bool{ return stk::simd::are_all(a); });
  this->set_expected_result([](){ return false; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, allOf_someTrue)
{
  if (this->is_scalar()) return;

  this->fill_alternating();
  this->compute_function([](const TypeParam& a)->bool{ return stk::simd::are_all(a); });
  this->set_expected_result([](){ return false; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, anyOf_allTrue)
{
  this->fill_constant(true);
  this->compute_function([](const TypeParam& a)->bool{ return stk::simd::are_any(a); });
  this->set_expected_result([](){ return true; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, anyOf_allFalse)
{
  this->fill_constant(false);
  this->compute_function([](const TypeParam& a)->bool{ return stk::simd::are_any(a); });
  this->set_expected_result([](){ return false; });
  this->verify();
}

TYPED_TEST(SimdBoolFixture, anyOf_someTrue)
{
  if (this->is_scalar()) return;

  this->fill_alternating();
  this->compute_function([](const TypeParam& a)->bool{ return stk::simd::are_any(a); });
  this->set_expected_result([](){ return true; });
  this->verify();
}
