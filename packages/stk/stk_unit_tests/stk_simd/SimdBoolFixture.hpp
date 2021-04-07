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

#ifndef SIMD_BOOL_FIXTURE_HPP 
#define SIMD_BOOL_FIXTURE_HPP 

#include <functional>

#include "stk_simd/Simd.hpp"

template<typename SimdType>
class SimdBoolFixture : public ::testing::Test
{
public:
  using FloatingPointType = typename stk::Traits<SimdType>::real_type;
  using FloatingPointScalar = typename stk::Traits<FloatingPointType>::base_type;

  bool is_scalar()
  {
    return SIZE == 1;
  }

  void fill_constant(bool a)
  {
    m_valueA = SimdType(a);
    
    load_scalar_values();
  }

  void fill_constant(bool a, bool b)
  {
    m_valueA = SimdType(a);
    m_valueB = SimdType(b);
    
    load_scalar_values();
  }

  void fill_alternating() 
  {
    FloatingPointType half(0.5);
    FloatingPointType zeroOne = zeros_ones();

    m_valueA = (zeroOne < half);
    m_valueB = (zeroOne < half);

    load_scalar_values();
  }

  void use_simd_value_as_result() 
  {
    m_result = m_valueA;
  }

  void set_expected_result(const std::function<bool()>& func)
  {
    for (int i=0; i<SIZE; i++) {
      m_expectedResult[i] = func();
    }
  }

  void set_expected_result(const std::function<bool(int)>& func)
  {
    for (int i=0; i<SIZE; i++) {
      m_expectedResult[i] = func(i);
    }
  }

  void compute_function(const std::function<SimdType(const SimdType&)>& func)
  {
    m_result = func(m_valueA);
  }

  void compute_function(const std::function<SimdType(const SimdType&, const SimdType&)>& func)
  {
    m_result = func(m_valueA, m_valueB);
  }

  void compute_expected_result(const std::function<bool(bool)>& func)
  {
    for (int i=0; i<SIZE; i++) {
      m_expectedResult[i] = func(m_valuesA[i]);
    }
  }

  void compute_expected_result(const std::function<bool(bool, bool)>& func)
  {
    for (int i=0; i<SIZE; i++) {
      m_expectedResult[i] = func(m_valuesA[i], m_valuesB[i]);
    }
  }

  void verify()
  {
    for (int i=0; i<SIZE; i++) {
      EXPECT_EQ(static_cast<bool>(m_result[i]), m_expectedResult[i]);
    }
  }

private:
  void load_scalar_values() 
  {
    for (int i=0; i<SIZE; i++) {
      m_valuesA[i] = static_cast<bool>(m_valueA[i]);
      m_valuesB[i] = static_cast<bool>(m_valueB[i]);
    }
  }

  FloatingPointType zeros_ones()
  {
    FloatingPointScalar data[SIZE];
    for (int i=0; i<SIZE; i++) {
      data[i] = i%2;
    }
    return stk::simd::load(data);
  }

  static constexpr int SIZE = stk::Traits<FloatingPointType>::length;

  bool m_valuesA[SIZE];
  bool m_valuesB[SIZE];

  SimdType m_valueA;
  SimdType m_valueB;

  bool m_expectedResult[SIZE];

  SimdType m_result;
};

#endif
