//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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
            setup_test<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::array_layout>(m, m / 10, repeat));
      },
      // run lambda will take the returned setup tuple
      [&](const int iteration, const int runsize, auto& data) {
        KokkosBlas::gemv("N", 1.0, data.A, data.x, 0.0, data.y);
      }));

  // return a vector of kernel base objects of type test_list
  return kernel_base_vector;
}
