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

#ifndef SIMD_LOAD_STORE_FIXTURE_HPP
#define SIMD_LOAD_STORE_FIXTURE_HPP

#include <gtest/gtest.h>
#include <array>

#include "stk_simd/Simd.hpp"

template <typename InputType>
class SimdLoadStoreFixture : public ::testing::Test
{
  public:
  using InputScalar = typename stk::Traits<InputType>::base_type;
  
  private:
  static constexpr int SIZE = stk::Traits<InputType>::length;
  static constexpr int N = 3;
  static constexpr int T_SIZE = N*SIZE;

  InputScalar data[T_SIZE];
 public:
  using SimdArray = std::array<InputType, N>;

  template <class Func>
  void fill(const Func& f)
  {
    for (int i = 0; i < T_SIZE; ++i) {
      data[i] = f(i);
    }
  }

  void fill_constant(const InputScalar& v)
  {
    fill([&v](int) { return v; });
  }

  template <class Func>
  void compare_until(int n, const Func& f) const
  {
    for (int i = 0; i < n; ++i) {
      EXPECT_EQ(f(i), data[i]) << i;
    }
  }

  template <class Func>
  void compare_from(int n, const Func& f) const
  {
    for (int i = n; i < T_SIZE; ++i) {
      EXPECT_EQ(f(i), data[i]) << i;
    }
  }

  template <class Func>
  void compare_until_with_offset(int n, const Func& f, const InputScalar& v) const
  {
    for(int i = 0; i < n; ++i) {
      int idx = i*N;
      EXPECT_EQ(f(idx), data[idx]) << idx;
      for(int k = 1; k < N; ++k) {
        EXPECT_EQ(v, data[idx+k]) << idx;
      }
    }
  }

  template <class Func>
  void compare_from_with_offset(int n, const Func& f) const
  {
    for (int i = n*N; i < T_SIZE; ++i) {
      EXPECT_EQ(f(i), data[i]) << i;
    }
  }

  template <class Func>
  void compare_until_array(int n, const Func& f) const
  {
    for (int i = 0; i < n*N; ++i) {
      EXPECT_EQ(f(i), data[i]) << i;
    }
  }

  template <class Func>
  void compare_from_array(int n, const Func& f) const
  {
    for (int i = n*N; i < T_SIZE; ++i) {
      EXPECT_EQ(f(i), data[i]) << i;
    }
  }

  InputType load_part(int numValid) const { return stk::simd::load_part(data, numValid); }
  InputType load_part_with_offset(int numValid) const { return stk::simd::load_part(data, N, numValid); }
  SimdArray load_array() const {
    SimdArray a;
    stk::simd::load_array<N>(a.data(), data);
    return a;
  }

  SimdArray load_array(int numValid) const {
    SimdArray a;
    stk::simd::load_array<N>(a.data(), data, numValid);
    return a;
  }

  void store_part(const InputType& src, int numValid) { stk::simd::store_part(data, src, numValid); }
  void store_part_with_offset(const InputType& src, int numValid) { stk::simd::store_part(data, src, N, numValid); }
  void store_array(const SimdArray& a) {
    stk::simd::store_array<N>(data, a.data());
  }
  void store_array(const SimdArray& a, int numValid) {
    stk::simd::store_array<N>(data, a.data(), numValid);
  }

  static constexpr int size() { return SIZE; }
  static constexpr int offset() {return N;}
};

#endif
