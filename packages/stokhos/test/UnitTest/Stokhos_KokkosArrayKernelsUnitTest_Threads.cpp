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

// Kernels
#include "Stokhos_Threads_CrsProductTensor.hpp"

// Tests
#include "Stokhos_KokkosArrayKernelsUnitTest.hpp"

using namespace KokkosKernelsUnitTest;

UnitTestSetup<Kokkos::Threads> setup;

// Test declarations
#include "Stokhos_KokkosArrayKernelsUnitTestDecl.hpp"
#include "Stokhos_KokkosArrayKernelsUnitTest_Host.hpp"

// Tests using Threads device
using Kokkos::Threads;
UNIT_TEST_GROUP_SCALAR_DEVICE( double, Threads )
UNIT_TEST_GROUP_SCALAR_HOST_DEVICE( double, Threads )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  const size_t team_count =
  Kokkos::hwloc::get_available_numa_count() *
    Kokkos::hwloc::get_available_cores_per_numa();
  const size_t threads_per_team =
    Kokkos::hwloc::get_available_threads_per_core();

  // Initialize threads
  Kokkos::InitializationSettings init_args;
  init_args.set_num_threads(team_count*threads_per_team);
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration( std::cout );

  // Setup (has to happen after initialization)
  setup.setup();

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
