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

#include "Kokkos_Macros.hpp"

#define SACADO_TEST_DFAD 1

#include "Fad_KokkosTests.hpp"

// Instantiate tests for Cuda device
using Kokkos::Cuda;
VIEW_FAD_TESTS_D( Cuda )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize Cuda
  Kokkos::finalize();

  return res;
}
