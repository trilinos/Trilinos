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

// Disable view specializations (changes atomic argument from ViewFadPtr to
// Fad*)
#define SACADO_DISABLE_FAD_VIEW_SPEC

#define SACADO_TEST_DFAD 1
#include "Fad_KokkosAtomicTests.hpp"

// Instantiate tests for OpenMP device

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
