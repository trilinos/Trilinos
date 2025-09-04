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

#include <gtest/gtest.h>

#include "SimdNameGenerator.hpp"
#include "SimdLoadStoreFixture.hpp"
#include "stk_simd/Simd.hpp"



template <class T>
using SimdLoadStore = SimdLoadStoreFixture<T>;
using SimdStoreTypes = ::testing::Types<stk::simd::Double, stk::simd::Float>;
TYPED_TEST_SUITE(SimdLoadStore, SimdStoreTypes, SimdNameGenerator);

TYPED_TEST(SimdLoadStore, partial_load_store)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;
  auto&& f = *this;

  for (int numValid = 0; numValid <= f.size(); ++numValid) {
    f.fill([](int i) { return static_cast<ScalarType>(i); });

    SimdType simd = f.load_part(numValid);
    simd *= 2.;

    f.fill_constant(-1.);
    f.store_part(simd, numValid);

    f.compare_until(numValid, [](int i) { return static_cast<ScalarType>(2 * i); });
    f.compare_from(numValid, [](int i) { return static_cast<ScalarType>(-1.); });
  }
}

TYPED_TEST(SimdLoadStore, partial_load_store_with_offset)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;
  auto&& f = *this;

  for (int numValid = 0; numValid <= f.size(); ++numValid) {
    f.fill([](int i) { return static_cast<ScalarType>(i); });

    SimdType simd = f.load_part_with_offset(numValid);
    simd *= 2.;

    f.fill_constant(-1.);
    f.store_part_with_offset(simd, numValid);

    f.compare_until_with_offset(numValid, [](int i) { return static_cast<ScalarType>(2 * i); }, -1.);
    f.compare_from_with_offset(numValid, [](int i) { return static_cast<ScalarType>(-1.); });
  }
}

TYPED_TEST(SimdLoadStore, array_load_store)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;
  auto&& f = *this;

  f.fill([](int i) { return static_cast<ScalarType>(i); });

  auto simdArray = f.load_array();
  for(auto&& simd : simdArray) {
    simd *= 2.;
  }

  f.fill_constant(-1.);
  f.store_array(simdArray);

  f.compare_from_array(0, [](int i) { return static_cast<ScalarType>(2*i); });
}

TYPED_TEST(SimdLoadStore, partial_array_load_store)
{
  using SimdType = TypeParam;
  using ScalarType = typename stk::Traits<SimdType>::base_type;
  auto&& f = *this;

  for (int numValid = 0; numValid <= f.size(); ++numValid) {
    f.fill([](int i) { return static_cast<ScalarType>(i); });

    auto simdArray = f.load_array(numValid);
    for(auto&& simd : simdArray) {
      simd *= 2.;
    }

    f.fill_constant(-1.);
    f.store_array(simdArray, numValid);

    f.compare_until_array(numValid, [](int i) { return static_cast<ScalarType>(2*i); });
    f.compare_from_array(numValid, [](int i) { return static_cast<ScalarType>(-1.); });
  }
}
