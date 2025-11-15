// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <iostream>
#include <cstdlib>
#include "Kokkos_Core.hpp"
#include "KokkosKernels_config.h"
#include "KokkosKernels_EagerInitialize.hpp"

int main() {
  Kokkos::initialize();
  // Schedule Kokkos::finalize to run when the program exits.
  std::atexit(Kokkos::finalize);
  // Now initialize all TPL singletons.
  // The test checks that program exits normally in this case.
  KokkosKernels::eager_initialize();
  return 0;
}
