// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#include <Kokkos_Core.hpp>
int main(int argc, char** argv)
{
  Kokkos::initialize(argc, argv);
  Kokkos::finalize();
  return 0;
}

