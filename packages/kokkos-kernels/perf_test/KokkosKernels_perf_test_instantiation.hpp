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
//
// Created by Berger-Vergiat, Luc on 2/6/23.
//

#ifndef KOKKOSKERNELS_PERF_TEST_INSTANTIATION_HPP
#define KOKKOSKERNELS_PERF_TEST_INSTANTIATION_HPP

#include "KokkosKernels_perf_test_utilities.hpp"

#ifndef KOKKOSKERNELS_PERF_TEST_NAME
#error "The macro KOKKOSKERNELS_PERF_TEST_NAME was not defined"
#endif

// All perf tests must implement print_options()
void print_options();

int main_instantiation(int argc, char** argv) {
  perf_test::CommonInputParams params;
  perf_test::parse_common_options(argc, argv, params);

  // If help is requested with "-h" or "--help", then just print the options
  // and quit.
  if (params.print_help) {
    print_options();
    return 0;
  }

  /* Assumption is that use_openmp/use_threads variables are */
  /* provided as numbers of threads */
  int num_threads = 1;
  if (params.use_openmp) {
    num_threads = params.use_openmp;
  } else if (params.use_threads) {
    num_threads = params.use_threads;
  }

  int device_id = 0;
  if (params.use_cuda)
    device_id = params.use_cuda - 1;
  else if (params.use_hip)
    device_id = params.use_hip - 1;
  else if (params.use_sycl)
    device_id = params.use_sycl - 1;

  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(num_threads).set_device_id(device_id));
  Kokkos::print_configuration(std::cout);
  std::cout << '\n';

  bool ran = false;

  if (params.use_openmp) {
#if defined(KOKKOS_ENABLE_OPENMP)
    std::cout << "Running on OpenMP backend.\n";
    KOKKOSKERNELS_PERF_TEST_NAME<Kokkos::OpenMP>(argc, argv, params);
    ran = true;
#else
    std::cout << "ERROR: OpenMP requested, but not available.\n";
    Kokkos::finalize();
    return 1;
#endif
  }
  if (params.use_threads) {
#if defined(KOKKOS_ENABLE_THREADS)
    std::cout << "Running on Threads backend.\n";
    KOKKOSKERNELS_PERF_TEST_NAME<Kokkos::Threads>(argc, argv, params);
    ran = true;
#else
    std::cout << "ERROR: Threads requested, but not available.\n";
    Kokkos::finalize();
    return 1;
#endif
  }
  if (params.use_cuda) {
#if defined(KOKKOS_ENABLE_CUDA)
    std::cout << "Running on Cuda backend.\n";
    KOKKOSKERNELS_PERF_TEST_NAME<Kokkos::Cuda>(argc, argv, params);
    ran = true;
#else
    std::cout << "ERROR: CUDA requested, but not available.\n";
    Kokkos::finalize();
    return 1;
#endif
  }
  if (params.use_hip) {
#if defined(KOKKOS_ENABLE_HIP)
    std::cout << "Running on HIP backend.\n";
    KOKKOSKERNELS_PERF_TEST_NAME<Kokkos::HIP>(argc, argv, params);
    ran = true;
#else
    std::cout << "ERROR: HIP requested, but not available.\n";
    Kokkos::finalize();
    return 1;
#endif
  }
  if (params.use_sycl) {
#if defined(KOKKOS_ENABLE_SYCL)
    std::cout << "Running on SYCL backend.\n";
    KOKKOSKERNELS_PERF_TEST_NAME<Kokkos::Experimental::SYCL>(argc, argv, params);
    ran = true;
#else
    std::cout << "ERROR: SYCL requested, but not available.\n";
    Kokkos::finalize();
    return 1;
#endif
  }
  if (!ran) {
#if defined(KOKKOS_ENABLE_SERIAL)
    std::cout << "Running on Serial backend.\n";
    KOKKOSKERNELS_PERF_TEST_NAME<Kokkos::Serial>(argc, argv, params);
#else
    std::cout << "ERROR: Tried to run on Serial device (as no parallel"
                 " backends requested), but Serial is not enabled.\n";
    Kokkos::finalize();
    return 1;
#endif
  }
  Kokkos::finalize();
  return 0;
}

// Undefine the macro to avoid potential bad interaction
// with other parts of the code...
#undef KOKKOSKERNELS_PERF_TEST_NAME

#endif  // KOKKOSKERNELS_PERF_TEST_INSTANTIATION_HPP
