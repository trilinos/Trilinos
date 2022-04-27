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
#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_Random.hpp>

// Required for tracked_testing version
#include "KokkosBlas3_gemm_tracked_perf_test.hpp"

#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>
#endif

// API ref for "General Matrix Multiplication" (gemm)
// https://github.com/kokkos/kokkos-kernels/wiki/BLAS-3%3A%3Agemm

// Usage: KokkosBlas::gemm(modeA, modeB, alpha, A, B, beta, C);
/*
* transA [in] "N" for non-transpose, "T" for transpose,
* "C" for conjugate transpose.
* All characters after the first are ignored. This works just like the BLAS
routines.
*
* transB [in] "N" for non-transpose,
* "T" for transpose,
* "C" for conjugate transpose.
* All characters after the first are ignored. This works just like the BLAS
routines.

* alpha [in] Input coefficient of A*x
* A [in] Input matrix, as a 2-D Kokkos::View
* B [in] Input matrix, as a 2-D Kokkos::View
* beta [in] Input coefficient of C
* C [in/out] Output vector, as a nonconst 2-D Kokkos::View

*/

// Define setup_test
template <class ExecSpace, class ALayout, class BLayout>
testData_gemm<ExecSpace, ALayout, BLayout> setup_test(int m, int n, int k,
                                                      int repeat) {
  testData_gemm<ExecSpace, ALayout, BLayout> testData_gemm_obj(m, n, k, repeat);

  return testData_gemm_obj;
}

test_list construct_gemm_kernel_base(const rajaperf::RunParams& run_params,
                                     const std::vector<gemmConfig>& m_n_k_vect)

{
  // instantiate kernel_base_vector as type test_list
  // kernel_base_vector will contain which tests to run, and data to run them
  test_list kernel_base_vector;

  for (const auto& value : m_n_k_vect) {
    kernel_base_vector.push_back(rajaperf::make_kernel_base(
        "BLAS3_GEMM_" + std::to_string(value.m) + "_" +
            std::to_string(value.n) + "_" + std::to_string(value.k),
        run_params,
        // setup_test lambda captures by value;
        // Mapping Kokkos features to RAJAPerf Suite
        // repeat = runreps (RAJAPerf Suite)
        // m = getActualRunSize() (RAJAPerf Suite)
        [=](const int repeat, const int) {
          // returns a tuple of testData objects
          return std::make_tuple(
              // setup_test is templated on ExecSpace and Layout
              setup_test<Kokkos::DefaultExecutionSpace,
                         Kokkos::DefaultExecutionSpace::array_layout,
                         Kokkos::DefaultExecutionSpace::array_layout>(
                  value.m, value.n, value.k, repeat));
        },
        // the run lambda will capture the returned setup_test tuple by
        // reference
        [&](const int iteration, const int runsize, auto& data) {
          // KokkosBlas::gemm(modeA, modeB, alpha, A, B, beta, C);
          KokkosBlas::gemm("N", "N", 1.0, data.A, data.B, 0.0, data.C);
        }));
  }

  // return a vector of kernel base objects of type test_list
  return kernel_base_vector;
}
