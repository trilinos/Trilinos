/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <KokkosBlas2_gemv.hpp>
#include <Kokkos_Random.hpp>

// Required for tracked_testing version
#include "KokkosBlas2_gemv_perf_test.hpp"
#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>
#endif

template <class ExecSpace, class Layout>
testData_gemv<ExecSpace, Layout> setup_test(int m, int n, int repeat) {
  testData_gemv<ExecSpace, Layout> testData_gemv_obj(m, n, repeat);

  return testData_gemv_obj;
}

test_list construct_gemv_kernel_base(const rajaperf::RunParams& run_params)

{
  // instantiate test_list as kernel_base_vector
  test_list kernel_base_vector;

  kernel_base_vector.push_back(rajaperf::make_kernel_base(
      "BLAS2_GEMV", run_params,
      // setup lambda by value
      [=](const int repeat, const int m) {
        // returns a tuple of testData_obj
        return std::make_tuple(
            setup_test<Kokkos::DefaultExecutionSpace,
                       Kokkos::DefaultExecutionSpace::array_layout>(m, m / 10,
                                                                    repeat));
      },
      // run lambda will take the returned setup tuple
      [&](const int iteration, const int runsize, auto& data) {
        KokkosBlas::gemv("N", 1.0, data.A, data.x, 0.0, data.y);
      }));

  // return a vector of kernel base objects of type test_list
  return kernel_base_vector;
}
