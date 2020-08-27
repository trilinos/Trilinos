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
#include <functional>

#include "stk_simd/Simd.hpp"


template<typename SimdType, typename ScalarType>
class SimdDoubleTwoArgs : public ::testing::Test
{
public:
  void fill_constant(double a, double b)
  {
    for (int i=0; i<stk::simd::ndoubles; i++) {
      m_valuesA[i] = a;
      m_valuesB[i] = b;
    }

    load_simd_values();
  }

  void fill_varying_linear_plus_minus()
  {
    for (int i=0; i<stk::simd::ndoubles; i++) {
      m_valuesA[i] = i;
      m_valuesB[i] = (i%2 == 0 ? 1.0 : -1.0);
    }

    load_simd_values();
  }

  void fill_varying_opposite_linear()
  {
    for (int i=0; i<stk::simd::ndoubles; i++) {
      m_valuesA[i] = i;
      m_valuesB[i] = stk::simd::ndoubles - i;
    }

    load_simd_values();
  }

  void compute_function(const std::function<SimdType(const stk::simd::Double&, const stk::simd::Double&)>& func)
  {
    m_result = func(m_valueA, m_valueB);
  }

  void compute_expected_result(const std::function<ScalarType(double, double)>& func)
  {
    for (int i=0; i<stk::simd::ndoubles; i++) {
      m_expectedResult[i] = func(m_valuesA[i], m_valuesB[i]);
    }
  }

  void verify()
  {
    for (int i=0; i<stk::simd::ndoubles; i++) {
      EXPECT_EQ(static_cast<ScalarType>(m_result[i]), m_expectedResult[i]);
    }
  }

private:
  void load_simd_values()
  {
    m_valueA = stk::simd::load(m_valuesA);
    m_valueB = stk::simd::load(m_valuesB);
  }

  double m_valuesA[stk::simd::ndoubles];
  double m_valuesB[stk::simd::ndoubles];

  stk::simd::Double m_valueA;
  stk::simd::Double m_valueB;

  ScalarType m_expectedResult[stk::simd::ndoubles];
  SimdType m_result;
};

using SimdDoubleMathTwoArgs = SimdDoubleTwoArgs<stk::simd::Double, double>;

TEST_F( SimdDoubleMathTwoArgs, copysign_posNeg )
{
  fill_constant(1.0, -2.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](double a, double b){ return std::copysign(a, b); });
  verify();
}

TEST_F( SimdDoubleMathTwoArgs, copysign_negPos )
{
  fill_constant(-2.0, 1.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](double a, double b){ return std::copysign(a, b); });
  verify();
}

TEST_F( SimdDoubleMathTwoArgs, copysign_negNeg )
{
  fill_constant(-3.0, -4.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](double a, double b){ return std::copysign(a, b); });
  verify();
}

TEST_F( SimdDoubleMathTwoArgs, copysign_posPos )
{
  fill_constant(6.0, 5.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](double a, double b){ return std::copysign(a, b); });
  verify();
}

TEST_F( SimdDoubleMathTwoArgs, copysign_varying )
{
  fill_varying_linear_plus_minus();
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return stk::math::copysign(a, b); });
  compute_expected_result([](double a, double b){ return std::copysign(a, b); });
  verify();
}

using SimdDoubleOperator = SimdDoubleTwoArgs<stk::simd::Bool, bool>;

TEST_F( SimdDoubleOperator, greaterThanEqual_aIsSmaller )
{
  fill_constant(1.0, 2.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a >= b; });
  compute_expected_result([](double a, double b){ return a >= b; });
  verify();
}

TEST_F( SimdDoubleOperator, greaterThanEqual_valuesAreEqual )
{
  fill_constant(2.0, 2.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a >= b; });
  compute_expected_result([](double a, double b){ return a >= b; });
  verify();
}

TEST_F( SimdDoubleOperator, greaterThanEqual_aIsLarger )
{
  fill_constant(3.0, -1.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a >= b; });
  compute_expected_result([](double a, double b){ return a >= b; });
  verify();
}

TEST_F( SimdDoubleOperator, greaterThanEqual_varying )
{
  fill_varying_opposite_linear();
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a >= b; });
  compute_expected_result([](double a, double b){ return a >= b; });
  verify();
}

TEST_F( SimdDoubleOperator, lessThan_aIsSmaller )
{
  fill_constant(1.0, 2.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a < b; });
  compute_expected_result([](double a, double b){ return a < b; });
  verify();
}

TEST_F( SimdDoubleOperator, lessThan_valuesAreEqual )
{
  fill_constant(2.0, 2.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a < b; });
  compute_expected_result([](double a, double b){ return a < b; });
  verify();
}

TEST_F( SimdDoubleOperator, lessThan_aIsLarger )
{
  fill_constant(3.0, -1.0);
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a < b; });
  compute_expected_result([](double a, double b){ return a < b; });
  verify();
}

TEST_F( SimdDoubleOperator, lessThan_varying )
{
  fill_varying_opposite_linear();
  compute_function([](const stk::simd::Double& a, const stk::simd::Double& b){ return a < b; });
  compute_expected_result([](double a, double b){ return a < b; });
  verify();
}

class SimdDoubleBool : public ::testing::Test
{
public:
  stk::simd::Double zeros_ones()
  {
    for (int i=0; i<stk::simd::ndoubles; i++) {
      m_data[i] = i%2;
    }
    return stk::simd::load(m_data);
  }

private:
  double m_data[stk::simd::ndoubles];
};

TEST_F(SimdDoubleBool, selectByLane_allFalse)
{
  stk::simd::Bool simdBool(false);

  for (int i=0; i<stk::simd::ndoubles; i++) {
    EXPECT_FALSE(simdBool[i]);
  }
}

TEST_F(SimdDoubleBool, selectByLane_allTrue)
{
  stk::simd::Bool simdTrue(true);

  for (int i=0; i<stk::simd::ndoubles; i++) {
    EXPECT_TRUE(simdTrue[i]);
  }
}

TEST_F(SimdDoubleBool, selectByLane_someTrue)
{
  stk::simd::Double half(0.5);
  stk::simd::Double zeroOne = zeros_ones();

  stk::simd::Bool simdBool = (zeroOne < half);

  for (int i=0; i<stk::simd::ndoubles; i++) {
    if (zeroOne[i] < half[i]) {
      EXPECT_TRUE(simdBool[i]);
    }
    else {
      EXPECT_FALSE(simdBool[i]);
    }
  }
}

TEST_F( SimdDoubleBool, operatorNot_allTrue )
{
  stk::simd::Bool value(true);
  stk::simd::Bool notValue = !value;

  for (int i=0; i<stk::simd::ndoubles; i++) {
    EXPECT_FALSE(notValue[i]);
  }
}

TEST_F( SimdDoubleBool, operatorNot_allFalse )
{
  stk::simd::Bool value(false);
  stk::simd::Bool notValue = !value;

  for (int i=0; i<stk::simd::ndoubles; i++) {
    EXPECT_TRUE(notValue[i]);
  }
}

TEST_F(SimdDoubleBool, operatorNot_someTrue)
{
  stk::simd::Double half(0.5);
  stk::simd::Double zeroOne = zeros_ones();

  stk::simd::Bool value = (zeroOne < half);
  stk::simd::Bool notValue = !value;

  for (int i=0; i<stk::simd::ndoubles; i++) {
    if (zeroOne[i] < half[i]) {
      EXPECT_FALSE(notValue[i]);
    }
    else {
      EXPECT_TRUE(notValue[i]);
    }
  }
}
