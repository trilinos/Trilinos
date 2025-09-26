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

#ifndef SIMD_FLOATING_POINT_FIXTURE_HPP 
#define SIMD_FLOATING_POINT_FIXTURE_HPP 

#include <functional>
#include <gtest/gtest.h>

#include "stk_simd/Simd.hpp"

template<typename InputType, typename ResultType>
class SimdFloatingPointFixture : public ::testing::Test
{
public:
  using InputScalar = typename stk::Traits<InputType>::base_type;
  using ResultScalar = typename stk::Traits<ResultType>::base_type;

  bool is_scalar()
  {
    return SIZE == 1;
  }

  void fill_constant(InputScalar a)
  {
    for (int i=0; i<SIZE; i++) {
      m_valuesA[i] = a;
    }

    load_simd_values();
  }

  void fill_constant(InputScalar a, InputScalar b)
  {
    for (int i=0; i<SIZE; i++) {
      m_valuesA[i] = a;
      m_valuesB[i] = b;
    }

    load_simd_values();
  }

  void fill_varying_linear_plus_minus()
  {
    for (int i=0; i<SIZE; i++) {
      m_valuesA[i] = i;
      m_valuesB[i] = (i%2 == 0 ? 1.0 : -1.0);
    }

    load_simd_values();
  }

  void fill_varying_opposite_linear()
  {
    for (int i=0; i<SIZE; i++) {
      m_valuesA[i] = i;
      m_valuesB[i] = SIZE - i;
    }

    load_simd_values();
  }

  void compute_function(const std::function<ResultType(const InputType&)>& func)
  {
    m_result = func(m_valueA);
  }

  void compute_function(const std::function<ResultType(const InputType&, const InputType&)>& func)
  {
    m_result = func(m_valueA, m_valueB);
  }

  void compute_expected_result(const std::function<ResultScalar(InputScalar)>& func)
  {
    for (int i=0; i<SIZE; i++) {
      m_expectedResult[i] = func(m_valuesA[i]);
    }
  }

  void compute_expected_result(const std::function<ResultScalar(InputScalar, InputScalar)>& func)
  {
    for (int i=0; i<SIZE; i++) {
      m_expectedResult[i] = func(m_valuesA[i], m_valuesB[i]);
    }
  }

  void verify()
  {
    for (int i=0; i<SIZE; i++) {
      EXPECT_EQ(static_cast<ResultScalar>(m_result[i]), m_expectedResult[i]);
    }
  }

private:
  void load_simd_values()
  {
    m_valueA = stk::simd::load(m_valuesA);
    m_valueB = stk::simd::load(m_valuesB);
  }

  static constexpr int SIZE = stk::Traits<InputType>::length;

  InputScalar m_valuesA[SIZE];
  InputScalar m_valuesB[SIZE];

  InputType m_valueA;
  InputType m_valueB;

  ResultScalar m_expectedResult[SIZE];

  ResultType m_result;
};

#endif
