// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_KokkosViewUQPCEUnitTest.hpp"

#include "Kokkos_Core.hpp"

// Instantiate test for Threads device
using Kokkos::Threads;
VIEW_UQ_PCE_TESTS_DEVICE( Threads )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize threads
  const size_t num_cores =
    Kokkos::hwloc::get_available_numa_count() *
    Kokkos::hwloc::get_available_cores_per_numa();
  const size_t num_hyper_threads =
    Kokkos::hwloc::get_available_threads_per_core();
  // const size_t num_cores = 1;
  // const size_t num_hyper_threads = 1;

  Kokkos::InitializationSettings init_args;
  init_args.set_num_threads(num_cores*num_hyper_threads);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
