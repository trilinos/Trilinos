// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file Test_Common_IndependentThreadScheduling.hpp
/// \brief Tests that KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS is set appropriately for known architectures

#ifndef TEST_COMMON_CUDAINDEPENDENTTHREADSCHEDULING
#define TEST_COMMON_CUDAINDEPENDENTTHREADSCHEDULING

#include <Kokkos_Core.hpp>
#include <KokkosKernels_Macros.hpp>

void test_cuda_independent_threads() {
#ifdef KOKKOS_ENABLE_CUDA

  // pull device off of default CUDA execution space
  const int device = Kokkos::Cuda{}.cuda_device();

  // retrieve device properties
  cudaDeviceProp deviceProp;
  cudaError_t err = cudaGetDeviceProperties(&deviceProp, device);
  ASSERT_EQ(err, cudaSuccess);  // no point in continuing if this fails

  // no independent thread scheduling for CC < 7.x
  if (deviceProp.major < 7) {
#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
    FAIL() << "KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS set, but found a GPU with compute capability < 70";
#endif
  } else if (7 == deviceProp.major) {  // VOLTA, TURING75
#ifndef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
    FAIL() << "KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS not set, but found a GPU with compute capability 7.x";
#endif
  } else if (8 == deviceProp.major) {  // AMPERE, ADA
#ifndef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
    FAIL() << "KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS not set, but found a GPU with compute capability 8.x";
#endif
  } else if (9 == deviceProp.major) {  // HOPPER
#ifndef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
    FAIL() << "KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS not set, but found a GPU with compute capability 9.x";
#endif
  } else if (10 == deviceProp.major && 0 == deviceProp.minor) {  // BLACKWELL
#ifndef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
    FAIL() << "KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS not set, but found a GPU with compute capability 10.x";
#endif
  } else {
    FAIL() << "Kokkos Kernels developer error, please report this. Ensure KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS is "
              "updated for this GPU architecture, and add the corresponding case to this test.";
  }

#else   // KOKKOS_ENABLE_CUDA
  GTEST_SKIP() << "CUDA backend not enabled.";
#endif  // KOKKOS_ENABLE_CUDA
}

TEST_F(TestCategory, common_cudaindependentthreads) { test_cuda_independent_threads(); }

#endif  // TEST_COMMON_CUDAINDEPENDENTTHREADSCHEDULING
