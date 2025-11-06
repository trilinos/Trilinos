// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef TEST_HIP_HPP
#define TEST_HIP_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_config.h>

#if defined(KOKKOSKERNELS_TEST_ETI_ONLY) && !defined(KOKKOSKERNELS_ETI_ONLY)
#define KOKKOSKERNELS_ETI_ONLY
#endif

class hip : public ::testing::Test {
 protected:
  static void SetUpTestCase() {}

  static void TearDownTestCase() {}
};

using HIPSpaceDevice        = Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>;
using HIPManagedSpaceDevice = Kokkos::Device<Kokkos::HIP, Kokkos::HIPManagedSpace>;

#define TestCategory hip

// Prefer <HIP, HIPSpace> for any testing where only one exec space is used
#if defined(KOKKOSKERNELS_INST_MEMSPACE_HIPMANAGEDSPACE) && !defined(KOKKOSKERNELS_INST_MEMSPACE_HIPSPACE)
#define TestDevice HIPManagedSpaceDevice
#else
#define TestDevice HIPSpaceDevice
#endif

#endif  // TEST_HIP_HPP
