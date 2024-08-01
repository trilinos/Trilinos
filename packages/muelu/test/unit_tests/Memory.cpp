// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Memory.hpp"

namespace {

TEUCHOS_UNIT_TEST(Memory, Basic) {
  using namespace MueLu::MemUtils;

  //    int n = 100000000;
  int n = 10000;

  out << "version: " << MueLu::Version() << std::endl;

  out << "Init:" << std::endl;

  out << " PrintMemoryUsage()" << PrintMemoryUsage() << std::endl;
  out << " PrintMemoryInfo()" << PrintMemoryInfo() << std::endl;

  out << "Allocation of " << n << " double (sizeof(double=" << sizeof(double) << "), total=" << sizeof(double) * n / (1024 * 1024) << " MB)" << std::endl;

  double* array = new double[n];

  // Use array to be sure that it will be allocated;
  for (int i = 0; i < n; i++)
    array[i] = i;

  out << " PrintMemoryUsage()" << PrintMemoryUsage() << std::endl;
  out << " PrintMemoryInfo()" << PrintMemoryInfo() << std::endl;

  out << "Delete:" << std::endl;

  delete[] array;

  out << " PrintMemoryUsage()" << PrintMemoryUsage() << std::endl;
  out << " PrintMemoryInfo()" << PrintMemoryInfo() << std::endl;

  success = true;
}

}  // namespace
