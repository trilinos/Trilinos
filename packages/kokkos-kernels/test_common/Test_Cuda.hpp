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
#ifndef TEST_CUDA_HPP
#define TEST_CUDA_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_config.h>

#if defined(KOKKOSKERNELS_TEST_ETI_ONLY) && !defined(KOKKOSKERNELS_ETI_ONLY)
#define KOKKOSKERNELS_ETI_ONLY
#endif

class Cuda : public ::testing::Test {
 protected:
  static void SetUpTestCase() {}

  static void TearDownTestCase() {}
};

#define TestCategory Cuda

using CudaSpaceDevice    = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
using CudaUVMSpaceDevice = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;

#ifdef KOKKOS_ENABLE_CUDA_UVM
// KOKKOS_ENABLE_CUDA_UVM macro and cmake option is deprecated
// But if it is defined, test with CudaUVMSpace.
// Make sure it's instantiated first:
#if defined(KOKKOSKERNELS_TEST_ETI_ONLY) && !defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE)
#error \
    "Deprecated option KOKKOS_ENABLE_CUDA_UVM is defined, so KokkosKernels will test with CudaUVMSpace. " \
  "KokkosKernels_INST_MEMSPACE_CUDAUVMSPACE=ON must be set in configuration."
#endif
#define TestDevice CudaUVMSpaceDevice
// Prefer <Cuda, CudaSpace> for any testing where only one exec space is used
#elif defined(KOKKOSKERNELS_INST_MEMSPACE_CUDAUVMSPACE) && !defined(KOKKOSKERNELS_INST_MEMSPACE_CUDASPACE)
#define TestDevice CudaUVMSpaceDevice
#else
#define TestDevice CudaSpaceDevice
#endif

#endif  // TEST_CUDA_HPP
