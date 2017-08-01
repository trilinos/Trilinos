// @HEADER
// ************************************************************************
//
//                           MiniTensor Package
//                 Copyright (2016) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#include <chrono>
#include <vector>

#include "MiniTensor.h"

struct Now {
  typedef std::chrono::time_point<std::chrono::high_resolution_clock> Impl;
  Impl impl;
};

Now now() {
  Now t;
  t.impl = std::chrono::high_resolution_clock::now();
  return t;
}

double operator-(Now b, Now a) {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(b.impl - a.impl)
             .count() *
         1e-9;
}

int
main(int ac, char* av[])
{
  Kokkos::initialize();

  int const n = 1000 * 1000;
  Kokkos::View<double*, Kokkos::DefaultExecutionSpace> results("results", n);
  auto t0 = now();
  for (int repeat_i = 0; repeat_i < 20; ++repeat_i) {
    Kokkos::parallel_for(Kokkos::RangePolicy<int, Kokkos::DefaultHostExecutionSpace>(0, n), [=](int i) {
      minitensor::Tensor<double, 3> B;
      B(0,0) = 3.14 * i;
      B(1,1) = 3.14 - double(i);
      B(2,2) = 3.14 + double(i);
      results(i) = trace(B);
    });
  }
  auto t1 = now();

  std::cout << "result(4242) " << results(4242) << '\n';
  std::cout << "total time " << (t1 - t0) << '\n';

  Kokkos::finalize();
}
