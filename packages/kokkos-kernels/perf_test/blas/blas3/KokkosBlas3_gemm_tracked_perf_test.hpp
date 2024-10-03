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

// This file is for the "tracked test" version of
// a Kokkos Kernels performance test.
// Created by David Poliakoff and Amy Powell on 9/22/2021

#ifndef KOKKOSKERNELS_KOKKOSBLAS_GEMM_TEST_RPS_HPP
#define KOKKOSKERNELS_KOKKOSBLAS_GEMM_TEST_RPS_HPP

#include <Kokkos_Core.hpp>
#include "blas/KokkosBlas3_gemm.hpp"
#include <Kokkos_Random.hpp>

// These headers are required for RPS tracked perf testing
#ifdef KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE
#include <PerfTestUtilities.hpp>

struct gemmConfig {
  int m;
  int n;
  int k;
};

test_list construct_gemm_kernel_base(const rajaperf::RunParams& run_params, const std::vector<gemmConfig>& n_k_vect);

#endif  // KOKKOSKERNELS_ENABLE_TESTS_AND_PERFSUITE

// Templating on these three types, mirroring
template <class ExecSpace, class ALayout, class BLayout>
struct testData_gemm {
  // Data for running tests
  // m is the number of rows in A
  int m = 1000;
  // n is the number of columns in A;
  // n is the number of rows in B;
  int n = 1000;
  // k is the number of columns in B;
  int k      = 1000;
  int repeat = 1;

  std::string modeA = "N";
  std::string modeB = "N";

  // Usage: KokkosBlas::gemm(modeA, modeB, alpha, A, B, beta, C);
  // Dense matrix-matrix multiply: C = beta*C + alpha*op(A)*op(B);
  // alpha:  alpha [in] Input coefficient of A*x
  // beta [in] Input coefficient of C

  double alpha = 1.0;
  double beta  = 0.0;

  using Scalar   = double;
  using MemSpace = typename ExecSpace::memory_space;
  using Device   = Kokkos::Device<ExecSpace, MemSpace>;

  // Create 2D Kokoks::View, "A" containing an input matrix with m x n
  // dimensions
  Kokkos::View<Scalar**, ALayout, Device> A;

  // Create 2D Kokkos::View, "B" containing an input matrix with n x k
  // dimensions
  Kokkos::View<Scalar**, BLayout, Device> B;

  // Create 2D Kokkos::View, "C" , the resultant matrix
  Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device> C;

  // class Constructor:
  testData_gemm(int m_in, int n_in, int k_in, int repeat_in) : m(m_in), n(n_in), k(k_in), repeat(repeat_in) {
    // You must set A, B and C equal to its intended value
    A = Kokkos::View<Scalar**, ALayout, Device>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), m, n);
    B = Kokkos::View<Scalar**, BLayout, Device>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "B"), n, k);
    C = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device>(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C"), m, k);
    // Seed random number generation
    Kokkos::Random_XorShift64_Pool<ExecSpace> pool(123);
    // Fill input matrices A and x with 10 random values from pool
    Kokkos::fill_random(A, pool, 10.0);
    Kokkos::fill_random(B, pool, 10.0);
  }
};

// Declare setup_test
template <typename ExecSpace, typename ALayout, typename BLayout>
testData_gemm<ExecSpace, ALayout, BLayout> setup_test(int m, int n, int k, int repeat);

#endif  // KOKKOSKERNELS_KOKKOSBLAS_GEMM_TEST_RPS_HPP
