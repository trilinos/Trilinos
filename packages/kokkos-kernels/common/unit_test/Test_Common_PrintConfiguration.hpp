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

/// \file Test_Common_PrintConfiguration.hpp
/// \brief Tests for print configuration

#ifndef KOKKOSKERNELS_PRINTCONFIGURATIONTEST_HPP
#define KOKKOSKERNELS_PRINTCONFIGURATIONTEST_HPP

#include "KokkosKernels_PrintConfiguration.hpp"

/// \brief Verify that all keys from kernels configuration and check their
/// values
void check_print_configuration(const std::ostringstream& msg) {
  bool kernelsVersionKeyFound   = false;
  bool enabledTPLsNamesKeyFound = false;
  // Iterate over lines returned from kokkos and extract key:value pairs
  std::stringstream ss{msg.str()};
  for (std::string line; std::getline(ss, line, '\n');) {
    auto found = line.find_first_of(':');
    if (found != std::string::npos) {
      auto currentKey = line.substr(0, found);
      if (currentKey == "  KokkosKernels Version") {
        kernelsVersionKeyFound = true;
      } else if (currentKey == "TPLs") {
        enabledTPLsNamesKeyFound = true;
      }
    }
  }
  EXPECT_TRUE(kernelsVersionKeyFound && enabledTPLsNamesKeyFound);
}

/// \brief Verify that print_configuration prints the expected keys from Kernels
/// configuration
template <typename exec_space>
void testPrintConfiguration() {
  // First, print this to cout in order to see what it looks like
  KokkosKernels::print_configuration(std::cout);
  // Then, run the actual test which prints the string to "out" and verifies
  // that out has meet some expected behavior
  std::ostringstream out;
  KokkosKernels::print_configuration(out);
  check_print_configuration(out);
}

TEST_F(TestCategory, common_print_configuration) { testPrintConfiguration<TestDevice>(); }

#endif  // KOKKOSKERNELS_PRINTCONFIGURATIONTEST_HPP
