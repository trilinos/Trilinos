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
