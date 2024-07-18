// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Utilities
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Device
#include "Kokkos_Core.hpp"

// Tests
#include "Stokhos_KokkosArrayKernelsUnitTest.hpp"

using namespace KokkosKernelsUnitTest;

UnitTestSetup<Kokkos::Serial> setup;

// Test declarations
#include "Stokhos_KokkosArrayKernelsUnitTestDecl.hpp"
#include "Stokhos_KokkosArrayKernelsUnitTest_Host.hpp"

// Tests using Serial device
using Kokkos::Serial;
UNIT_TEST_GROUP_SCALAR_DEVICE( double, Serial )
UNIT_TEST_GROUP_SCALAR_HOST_DEVICE( double, Serial )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize serial
  Kokkos::initialize();

  // Setup (has to happen after initialization)
  setup.setup();

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
