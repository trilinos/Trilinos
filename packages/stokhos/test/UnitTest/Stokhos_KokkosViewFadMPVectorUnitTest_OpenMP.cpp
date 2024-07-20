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

#include "Stokhos_KokkosViewFadMPVectorUnitTest.hpp"

#include "Kokkos_Core.hpp"

// Instantiate test for OpenMP device
using Kokkos::OpenMP;
VIEW_FAD_MP_VECTOR_TESTS_DEVICE( OpenMP )

int main( int argc, char* argv[] ) {
// Setup the MPI session
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Kokkos
  Kokkos::initialize(argc, argv);

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
