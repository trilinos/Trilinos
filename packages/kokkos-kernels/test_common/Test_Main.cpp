// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);

  ::testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  Kokkos::finalize();

  return result;
}
