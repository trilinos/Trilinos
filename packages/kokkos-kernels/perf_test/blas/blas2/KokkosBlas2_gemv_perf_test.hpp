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

#ifndef KOKKOSKERNELS_KOKKOSBLAS_GEMV_TEST_RPS_HPP
#define KOKKOSKERNELS_KOKKOSBLAS_GEMV_TEST_RPS_HPP

#include <Kokkos_Core.hpp>
#include "blas/KokkosBlas1_dot.hpp"
#include <Kokkos_Random.hpp>

// These headers are required for RPS perf test implementation
//
#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>

test_list construct_gemv_kernel_base(const rajaperf::RunParams& run_params);

#endif  // KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE

template <class ExecSpace, class Layout>
struct testData_gemv {
  // type aliases for Kokkos data structures

  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  // These are fields in the struct
  // m is vector length
  int m = 100000;
  // n is the number of columns
  int n           = 100000;
  int repeat      = 1;
  bool layoutLeft = true;

  // Create 2D view, "A", w/ Device as the ExecSpace; this is the input matrix,
  // as a 2-D Kokkos::View
  Kokkos::View<Scalar**, Layout, Device> A;

  // Create 1D view, "x", w/ Device as the ExecSpace; input vector, as a 1-D
  // Kokkos::View
  Kokkos::View<Scalar*, Device> x;

  // Create 1D view, "y", w/ Device as the ExecSpace; input/output vector, as a
  // 1-D Kokkos::View
  Kokkos::View<Scalar*, Device> y;

  // A function with no return type whose name is the name of the class is a
  // constructor or a destructor;
  // Constructor -- create function:
  testData_gemv(int m_in, int n_in, int repeat_in) : m(m_in), n(n_in), repeat(repeat_in) {
    A = Kokkos::View<Scalar**, Layout, Device>(Kokkos::ViewAllocateWithoutInitializing("A"), m, n);
    x = Kokkos::View<Scalar*, Device>(Kokkos::ViewAllocateWithoutInitializing("x"), n);
    y = Kokkos::View<Scalar*, Device>(Kokkos::ViewAllocateWithoutInitializing("y"), m);

    Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);

    // Fill input matrices A and x with 10 random values from pool
    Kokkos::fill_random(A, pool, 10.0);
    Kokkos::fill_random(x, pool, 10.0);
  }
};

template <typename ExecSpace, typename Layout>
testData_gemv<ExecSpace, Layout> setup_test(int m, int n, int repeat, bool layoutLeft);

#endif  // KOKKOSKERNELS_KOKKOSBLAS_GEMV_TEST_RPS_HPP
