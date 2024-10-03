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

// Created by David Poliakoff and Amy Powell on 6/15/2021

#ifndef KOKKOSKERNELS_KOKKOSBLAS_DOT_TEST_RPS_HPP
#define KOKKOSKERNELS_KOKKOSBLAS_DOT_TEST_RPS_HPP

#include <Kokkos_Core.hpp>
#include "KokkosBlas1_dot.hpp"
#include <Kokkos_Random.hpp>

// These headers are required for RPS perf test implementation
//
#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>
test_list construct_dot_kernel_base(const rajaperf::RunParams& run_params);
#endif  // KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE

template <class ExecSpace>
struct testData {
  // type aliases
  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  // Run Time Info for KK implementation
  //  int use_cuda    = 0;
  //  int use_openmp  = 0;
  //  int use_threads = 0;

  // m is vector length
  int m      = 100000;
  int repeat = 1;

  // Test Matrices x and y, View declaration

  // Create 1D view w/ Device as the ExecSpace; this is an input vector
  Kokkos::View<Scalar*, Device> x;

  // Create 1D view w/ Device as the ExecSpace; this is the output vector
  Kokkos::View<Scalar*, Device> y;

  // A function with no return type whose name is the name of the class is a
  // constructor or a destructor;
  // Constructor -- create function:
  testData(int m_in) : m(m_in) {
    x = Kokkos::View<Scalar*, Device>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "x"), m);
    y = Kokkos::View<Scalar*, Device>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "y"), m);

    Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);

    Kokkos::fill_random(x, pool, 10.0);
    Kokkos::fill_random(y, pool, 10.0);
  }
};

template <typename ExecSpace>
testData<ExecSpace> setup_test(int m, int repeat);

#endif  // KOKKOSKERNELS_KOKKOSBLAS_DOT_TEST_HPP
