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
#include <KokkosBlas1_dot.hpp>
#include <Kokkos_Random.hpp>

// For RPS implementation
#include "KokkosBlas_dot_perf_test.hpp"
#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>
#endif

template <class ExecSpace>
testData<ExecSpace> setup_test(int m, int repeat) {
  // use constructor to generate test data
  testData<ExecSpace> testData_obj(m);

  // set a field in the struct
  testData_obj.m      = m;
  testData_obj.repeat = repeat;

  return testData_obj;
}

test_list construct_dot_kernel_base(const rajaperf::RunParams& run_params)

{
  // instantiate test_list as kernel_base_vector
  test_list kernel_base_vector;

  kernel_base_vector.push_back(rajaperf::make_kernel_base(
      "BLAS_DOT ", run_params,
      [=](const int repeat, const int m) {
        // returns a tuple of testData_obj
        return std::make_tuple(setup_test<Kokkos::DefaultExecutionSpace>(m, repeat));
      },
      [&](const int, const int, auto& data) { KokkosBlas::dot(data.x, data.y); }));

  // return a vector of kernel base objects of type test_list
  return kernel_base_vector;
}
