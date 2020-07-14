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

#include "stk_unit_test_utils/PrintType.hpp"

#ifndef STK_KOKKOS_SIMD
#define STK_KOKKOS_SIMD
#endif

#ifndef USE_STK_SIMD_NONE
#define USE_STK_SIMD_NONE
#endif

#include "stk_simd/Simd.hpp"

TEST( SimdInfoScalar, printWidths )
{
  std::cout << "stk::simd::nfloats " << stk::simd::nfloats << std::endl;
  std::cout << "stk::simd::ndoubles " << stk::simd::ndoubles << std::endl;
}

TEST( SimdInfoScalar, printTypes )
{
  std::cout << "Datatype stored by stk::simd::Float is ";
  stk::simd::Float f;
  stk::unit_test_util::print_type(f._data);

  std::cout << "Datatype stored by stk::simd::Double is ";
  stk::simd::Double d;
  stk::unit_test_util::print_type(d._data);
}

using FloatDataScalar = SIMD_NAMESPACE::simd<float, SIMD_NAMESPACE::simd_abi::scalar>;
using DoubleDataScalar = SIMD_NAMESPACE::simd<double, SIMD_NAMESPACE::simd_abi::scalar>;

TEST( SimdInfoScalar, checkTypes )
{
  stk::simd::Float f;
  EXPECT_TRUE((std::is_same<decltype(f._data), FloatDataScalar>::value));

  stk::simd::Double d;
  EXPECT_TRUE((std::is_same<decltype(d._data), DoubleDataScalar>::value));
}
