// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
// 


#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t
#include <stk_simd/Simd.hpp>
#include <limits>                       // for std::numeric_limits

TEST(Simd, basic)
{
#if defined(STK_SIMD_AVX512) || defined(STK_SIMD_AVX) || defined(STK_SIMD_SSE)
    EXPECT_EQ(stk::simd::nfloats, 2*stk::simd::ndoubles);
#endif
}

TEST(Simd, whichInstructions)
{
#if defined(STK_SIMD_AVX512)
   std::cout<<"STK_SIMD_AVX512";
#elif defined(STK_SIMD_AVX)
   std::cout<<"STK_SIMD_AVX";
#elif defined(STK_SIMD_SSE)
   std::cout<<"STK_SIMD_SSE";
#else
   std::cout<<"no simd instructions!"<<std::endl;
#endif
   std::cout<<", stk::simd::ndoubles="<<stk::simd::ndoubles<<std::endl;
}

//BEGINsimdSimd
TEST(stkMeshHowTo, simdSimdTest)
{
  const int N = 512; // this is a multiple of the simd width
                     // if this is not true, the remainder 
                     // must be handled appropriately
                     
  static_assert( N % stk::simd::ndoubles == 0, "Required to be a multiple of ndoubles");

  std::vector<double, non_std::AlignedAllocator<double,64> > x(N);
  std::vector<double, non_std::AlignedAllocator<double,64> > y(N);
  std::vector<double, non_std::AlignedAllocator<double,64> > solution(N);
  
  for (int n=0; n < N; ++n) {
    x[n] = (rand()-0.5)/RAND_MAX;
    y[n] = (rand()-0.5)/RAND_MAX;
  }

  for (int n=0; n < N; n+=stk::simd::ndoubles) {
    const stk::simd::Double xl = stk::simd::load(&x[n]);
    const stk::simd::Double yl = stk::simd::load(&y[n]);
    stk::simd::Double zl = stk::math::abs(xl) * stk::math::exp(yl);
    stk::simd::store(&solution[n],zl);
  }  

  const double epsilon = std::numeric_limits<double>::epsilon();
  for (int n=0; n < N; ++n) {
    EXPECT_NEAR( std::abs(x[n]) * std::exp(y[n]), solution[n], epsilon );
  }
}
//ENDsimdSimd
