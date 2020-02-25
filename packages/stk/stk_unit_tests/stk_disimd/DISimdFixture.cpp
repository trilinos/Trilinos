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

#include "DISimdFixture.hpp"

namespace stk
{
namespace unit_test_simd
{

std::string get_simd_abi_native_name()
{
#if defined(__CUDACC__)
  return "scalar (cuda)";
#elif defined(__AVX512F__)
  return "avx512";
#elif defined(__AVX__)
  return "avx";
#elif defined(__SSE2__)
  return "sse";
#elif defined(__ARM_NEON)
  return "neon (ARM)";
#elif defined(__VSX__)
  return "vsx (IBM POWER)";
#else
  return "pack<16> (no simd defined)";
#endif
}

std::vector<double> X_RandomValues(int N)
{
  std::vector<double> x(N);
  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.4)/RAND_MAX;
  }
  return x;
}

std::vector<double> Y_RandomValues(std::vector<double>& x)
{
  int N = x.size();
  std::vector<double> y(N);
  for (int n=0; n < N; ++n) {
    y[n] = 0.5*x[n];
  }
  return y;
}

}
}
