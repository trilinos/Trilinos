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

#ifndef STK_KOKKOS_SIMD
#define STK_KOKKOS_SIMD
#endif

#include "Kokkos_Core.hpp"
#include "stk_simd/Simd.hpp"

int get_float_width_on_device()
{
  int result;

  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const int& i, int& width) {
    stk::simd::Float f;
    width = f._data.size();
  }, result);

  return result;
}

int get_double_width_on_device()
{
  int result;

  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const int& i, int& width) {
    stk::simd::Double d;
    width = d._data.size();
  }, result);

  return result;
}

#ifdef KOKKOS_ENABLE_CUDA
TEST( SimdInfoDevice, checkWidths )
{
  EXPECT_EQ(get_float_width_on_device(), 1);
  EXPECT_EQ(get_double_width_on_device(), 1);
}

TEST( SimdInfoDevice, printWidths )
{
  std::cout << "width of stk::simd::Float on device " << get_float_width_on_device() << std::endl;
  std::cout << "width of stk::simd::Double on device " << get_double_width_on_device() << std::endl;
}
#endif
