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

using SimdFloatBoolOperator = SimdBoolFixture<stk::simd::Boolf, stk::simd::nfloats>;

TEST_F(SimdFloatBoolOperator, selectByLane_allTrue)
{
  fill_constant(true);
  use_simd_value_as_result();
  set_expected_result([](int i){ return true; });
  verify();
}

TEST_F(SimdFloatBoolOperator, selectByLane_allFalse)
{
  fill_constant(false);
  use_simd_value_as_result();
  set_expected_result([](int i){ return false; });
  verify();
}

TEST_F(SimdFloatBoolOperator, selectByLane_someTrue)
{
  if (is_scalar()) return;

  fill_alternating();
  use_simd_value_as_result();
  set_expected_result([](int i){ return (i%2 == 0); });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorOr_TrueTrue)
{
  fill_constant(true, true);
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a || b; });
  compute_expected_result([](bool a, bool b){ return a || b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorOr_FalseFalse)
{
  fill_constant(false, false);
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a || b; });
  compute_expected_result([](bool a, bool b){ return a || b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorOr_TrueFalse)
{
  fill_constant(true, false);
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a || b; });
  compute_expected_result([](bool a, bool b){ return a || b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorOr_Alternating)
{
  if (is_scalar()) return;

  fill_alternating();
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a || b; });
  compute_expected_result([](bool a, bool b){ return a || b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorAnd_TrueTrue)
{
  fill_constant(true, true);
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a && b; });
  compute_expected_result([](bool a, bool b){ return a && b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorAnd_FalseFalse)
{
  fill_constant(false, false);
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a && b; });
  compute_expected_result([](bool a, bool b){ return a && b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorAnd_TrueFalse)
{
  fill_constant(true, false);
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a && b; });
  compute_expected_result([](bool a, bool b){ return a && b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorAnd_Alternating)
{
  if (is_scalar()) return;

  fill_alternating();
  compute_function([](const stk::simd::Boolf& a, const stk::simd::Boolf& b){ return a && b; });
  compute_expected_result([](bool a, bool b){ return a && b; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorNot_allTrue)
{
  fill_constant(true);
  compute_function([](const stk::simd::Boolf& a){ return !a; });
  compute_expected_result([](bool a){ return !a; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorNot_allFalse)
{
  fill_constant(false);
  compute_function([](const stk::simd::Boolf& a){ return !a; });
  compute_expected_result([](bool a){ return !a; });
  verify();
}

TEST_F(SimdFloatBoolOperator, operatorNot_someTrue)
{
  if (is_scalar()) return;

  fill_alternating();
  compute_function([](const stk::simd::Boolf& a){ return !a; });
  compute_expected_result([](bool a){ return !a; });
  verify();
}

TEST_F(SimdFloatBoolOperator, allOf_allTrue)
{
  fill_constant(true);
  compute_function([](const stk::simd::Boolf& a)->bool{ return stk::simd::are_all(a); });
  set_expected_result([](){ return true; });
  verify();
}

TEST_F(SimdFloatBoolOperator, allOf_allFalse)
{
  fill_constant(false);
  compute_function([](const stk::simd::Boolf& a)->bool{ return stk::simd::are_all(a); });
  set_expected_result([](){ return false; });
  verify();
}

TEST_F(SimdFloatBoolOperator, allOf_someTrue)
{
  if (is_scalar()) return;

  fill_alternating();
  compute_function([](const stk::simd::Boolf& a)->bool{ return stk::simd::are_all(a); });
  set_expected_result([](){ return false; });
  verify();
}

TEST_F(SimdFloatBoolOperator, anyOf_allTrue)
{
  fill_constant(true);
  compute_function([](const stk::simd::Boolf& a)->bool{ return stk::simd::are_any(a); });
  set_expected_result([](){ return true; });
  verify();
}

TEST_F(SimdFloatBoolOperator, anyOf_allFalse)
{
  fill_constant(false);
  compute_function([](const stk::simd::Boolf& a)->bool{ return stk::simd::are_any(a); });
  set_expected_result([](){ return false; });
  verify();
}

TEST_F(SimdFloatBoolOperator, anyOf_someTrue)
{
  if (is_scalar()) return;

  fill_alternating();
  compute_function([](const stk::simd::Boolf& a)->bool{ return stk::simd::are_any(a); });
  set_expected_result([](){ return true; });
  verify();
}
