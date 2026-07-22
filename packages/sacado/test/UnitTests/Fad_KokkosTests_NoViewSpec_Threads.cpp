// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Disable view specializations
#define SACADO_DISABLE_FAD_VIEW_SPEC

#define SACADO_TEST_DFAD 1
#include "Fad_KokkosTests.hpp"

// Instantiate tests for Threads device
using Kokkos::Threads;
VIEW_FAD_TESTS_D( Threads )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize threads
  size_t num_cores =
    Kokkos::hwloc::get_available_numa_count() *
    Kokkos::hwloc::get_available_cores_per_numa();
  size_t num_hyper_threads =
    Kokkos::hwloc::get_available_threads_per_core();

  Kokkos::InitializationSettings init_args;
  init_args.set_num_threads(num_cores*num_hyper_threads);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize threads
  Kokkos::finalize();

  return res;
}
