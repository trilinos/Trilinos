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

#include "Fad_Fad_KokkosTests.hpp"

#include "Kokkos_Core.hpp"

// Instantiate tests for OpenMP device
#define SACADO_TEST_DFAD 1
using Kokkos::OpenMP;
VIEW_FAD_TESTS_D( OpenMP )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize OpenMP
  Kokkos::initialize();
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize OpenMP
  Kokkos::finalize();

  return res;
}
