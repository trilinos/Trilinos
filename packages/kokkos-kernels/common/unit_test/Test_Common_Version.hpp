// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file Test_Common_Version.hpp
/// \brief Tests that the version information that Kokkos Kernels
///        makes available in KokkosKernels_config.h is properly
///        accessible and correct.

#ifndef TEST_COMMON_VERSION_HPP
#define TEST_COMMON_VERSION_HPP

#include <Kokkos_Core.hpp>
#include <KokkosKernels_config.h>

void test_version_info() {
#ifndef KOKKOSKERNELS_VERSION
  static_assert(false, "KOKKOSKERNELS_VERSION macro is not defined!");
#endif

#ifndef KOKKOSKERNELS_VERSION_MAJOR
  static_assert(false, "KOKKOSKERNELS_VERSION_MAJOR macro is not defined!");
#endif

#ifndef KOKKOSKERNELS_VERSION_MINOR
  static_assert(false, "KOKKOSKERNELS_VERSION_MINOR macro is not defined!");
#endif

#ifndef KOKKOSKERNELS_VERSION_PATCH
  static_assert(false, "KOKKOSKERNELS_VERSION_PATCH macro is not defined!");
#endif

  static_assert(KOKKOSKERNELS_VERSION == (KOKKOSKERNELS_VERSION_MAJOR * 10000 + KOKKOSKERNELS_VERSION_MINOR * 100 +
                                          KOKKOSKERNELS_VERSION_PATCH));
}

TEST_F(TestCategory, common_version) { test_version_info(); }

#endif  // TEST_COMMON_VERSION_HPP
