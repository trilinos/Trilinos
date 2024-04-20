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

#include "SimdDeviceWidths.hpp"
#include <stk_simd_view/simd_parallel.hpp>

using ExeHost = Kokkos::DefaultHostExecutionSpace;
using ExeDev = Kokkos::DefaultExecutionSpace;

template <typename ExecSpace>
struct FuncWithExecPolicy {
  using execution_policy = Kokkos::RangePolicy<ExecSpace>;
};

#ifdef STK_ENABLE_GPU
TEST( SimdInfo, checkDeviceWidths )
{
  EXPECT_EQ(stk::unit_test_util::get_float_width_on_device(), 1);
  EXPECT_EQ(stk::unit_test_util::get_double_width_on_device(), 1);
}

void do_device_test()
{
  const int N = 5;
  int returnedVal = 0;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int i, int& localVal)
  {
    localVal = stk::simd::get_simd_loop_size<float,FuncWithExecPolicy<ExeDev>>(N);
  }, returnedVal);
  EXPECT_EQ(N, returnedVal);
}

TEST( SimdInfo, checkDeviceLoopLength )
{
  do_device_test();
}

#endif

using FloatDataNative = SIMD_NAMESPACE::simd<float, SIMD_NAMESPACE::simd_abi::native>;
using DoubleDataNative = SIMD_NAMESPACE::simd<double, SIMD_NAMESPACE::simd_abi::native>;

TEST( SimdInfo, checkTypes )
{
  stk::simd::Float f;
  EXPECT_TRUE((std::is_same<decltype(f._data), FloatDataNative>::value));

  stk::simd::Double d;
  EXPECT_TRUE((std::is_same<decltype(d._data), DoubleDataNative>::value));
}

TEST( SimdInfo, checkHostLoopLength )
{
  const int N = 19;
  int loopSize = stk::simd::get_simd_loop_size<float,FuncWithExecPolicy<ExeHost>>(N);
  const int expectedLoopSize = (N%stk::simd::nfloats==0 ? 0 : 1) + N/stk::simd::nfloats;
  EXPECT_EQ(expectedLoopSize, loopSize);
}

