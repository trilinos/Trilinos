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

#include <functional>
#include <type_traits>
#include <gtest/gtest.h>

#include "stk_simd/Simd.hpp"

#include <cstdlib>

template<typename InputType, typename ResultType, typename InputTypeB = InputType>
class SimdFloatingPointFixture : public ::testing::Test
{
public:
  using InputScalarA = typename stk::Traits<InputType>::base_type;
  using InputScalarB = typename stk::Traits<InputTypeB>::base_type;
  using ResultScalar = typename stk::Traits<ResultType>::base_type;

  static constexpr int SIZE_A = stk::Traits<InputType>::length;
  static constexpr int SIZE_B = stk::Traits<InputTypeB>::length;
  static constexpr int SIZE_R = stk::Traits<ResultType>::length;

  static_assert(SIZE_B == 1 || SIZE_B == SIZE_A,
    "InputTypeB must be scalar (length 1) or same SIMD width as InputType.");

  bool is_scalar() const { return SIZE_A == 1; }
  bool is_scalar_a() const { return SIZE_A == 1; }
  bool is_scalar_b() const { return SIZE_B == 1; }

  // Fill A with a constant (existing behavior preserved).
  void fill_constant(InputScalarA a)
  {
    for (int i = 0; i < SIZE_A; ++i) m_valuesA[i] = a;
    load_simd_values();
  }

  void fill_constant(InputScalarA a, InputScalarB b)
  {
    for (int i = 0; i < SIZE_A; ++i) m_valuesA[i] = a;
    if constexpr (SIZE_B == 1) {
      m_valuesB[0] = b;
    } else {
      for (int i = 0; i < SIZE_B; ++i) m_valuesB[i] = b;
    }

    load_simd_values();
  }

  void fill_constant_b(InputScalarB b)
  {
    if constexpr (SIZE_B == 1) {
      m_valuesB[0] = b;
    } else {
      for (int i = 0; i < SIZE_B; ++i) m_valuesB[i] = b;
    }
    load_simd_values();
  }

  void fill_varying_linear_plus_minus()
  {
    for (int i = 0; i < SIZE_A; ++i) m_valuesA[i] = static_cast<InputScalarA>(i);
    if constexpr (SIZE_B == 1) {
      m_valuesB[0] = static_cast<InputScalarB>(1.0);
    } else {
      for (int i = 0; i < SIZE_B; ++i)
        m_valuesB[i] = static_cast<InputScalarB>((i % 2 == 0) ? 1.0 : -1.0);
    }

    load_simd_values();
  }

  void fill_varying_opposite_linear()
  {
    for (int i = 0; i < SIZE_A; ++i) m_valuesA[i] = static_cast<InputScalarA>(i);
    if constexpr (SIZE_B == 1) {
      m_valuesB[0] = static_cast<InputScalarB>(SIZE_A);
    } else {
      for (int i = 0; i < SIZE_B; ++i) m_valuesB[i] = static_cast<InputScalarB>(SIZE_A - i);
    }

    load_simd_values();
  }

  void fill_linear()
  {
    // so that we don't just test 0,1 for SSE.
    const int initial_a = 2;
    const int initial_b = 3;

    for (int i = 0; i < SIZE_A; ++i) {
      m_valuesA[i] = static_cast<InputScalarA>(initial_a + i);
    }
    if constexpr (SIZE_B == 1) {
      m_valuesB[0] = InputScalarB(42);
    } else {
      for (int i = 0; i < SIZE_B; ++i) {
        m_valuesA[i] = static_cast<InputScalarB>(initial_b + i);
      }
    }

    load_simd_values();
  }

  void compute_function(const std::function<ResultType(const InputType&)>& func)
  {
    m_result = func(m_valueA);
  }

  void compute_function(const std::function<ResultType(const InputType&, const InputTypeB&)>& func)
  {
    m_result = func(m_valueA, m_valueB);
  }

  void compute_expected_result(const std::function<ResultScalar(InputScalarA)>& func)
  {
    for (int i = 0; i < SIZE_A; ++i) {
      m_expectedResult[i] = func(m_valuesA[i]);
    }
  }

  void compute_expected_result(const std::function<ResultScalar(InputScalarA, InputScalarB)>& func)
  {
    for (int i = 0; i < SIZE_A; ++i) {
      const InputScalarA a = m_valuesA[i];
      const InputScalarB b = (SIZE_B == 1) ? m_valuesB[0] : m_valuesB[i];
      m_expectedResult[i] = func(a, b);
    }
  }

  void verify()
  {
    for (int i = 0; i < SIZE_A; ++i) {
      // lane-wise compare: ResultType is expected to be SIMD with SIZE_A lanes
      if(std::is_same_v<ResultScalar,float>) {
        EXPECT_FLOAT_EQ(static_cast<ResultScalar>(m_result[i]), m_expectedResult[i]);
      } else if (std::is_same_v<ResultScalar,double>) {
        EXPECT_DOUBLE_EQ(static_cast<ResultScalar>(m_result[i]), m_expectedResult[i]);
      } else {
        EXPECT_EQ(static_cast<ResultScalar>(m_result[i]), m_expectedResult[i]);
      }
    }
  }

  void verify_with_tolerance(const double tolerance = 1e-12)
  {
    for (int i = 0; i < SIZE_A; ++i) {
      // lane-wise compare: ResultType is expected to be SIMD with SIZE_A lanes
      ASSERT_NEAR(static_cast<ResultScalar>(m_result[i]), m_expectedResult[i], tolerance);
    }
  }

private:
  // load A and B with proper handling for scalar-vs-SIMD B
  void load_simd_values()
  {
    m_valueA = stk::simd::load(m_valuesA);

    if constexpr (SIZE_B == 1) {
      // For scalar B (e.g., int/double), just store the single value; no SIMD load.
      m_valueB = m_valuesB[0];
    } else {
      m_valueB = stk::simd::load(m_valuesB);
    }
  }

  // Storage (A and B can have different scalar types/widths)
  InputScalarA m_valuesA[SIZE_A] = {};
  InputScalarB m_valuesB[(SIZE_B == 0 ? 1 : SIZE_B)] = {}; // guard for weird traits

  InputType   m_valueA{};
  InputTypeB  m_valueB{};

  ResultScalar m_expectedResult[SIZE_A] = {};
  ResultType   m_result{};
};

