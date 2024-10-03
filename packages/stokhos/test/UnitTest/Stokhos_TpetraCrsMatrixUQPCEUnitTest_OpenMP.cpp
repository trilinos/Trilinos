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

#include "Stokhos_TpetraCrsMatrixUQPCEUnitTest.hpp"

#include "Kokkos_Core.hpp"
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

// Instantiate tests for threads node
typedef Tpetra::KokkosCompat::KokkosDeviceWrapperNode<Kokkos::OpenMP> OpenMPWrapperNode;
CRSMATRIX_UQ_PCE_TESTS_N( OpenMPWrapperNode )

int main( int argc, char* argv[] ) {
// Setup the MPI session
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Kokkos
  Kokkos::initialize(argc, argv);

  // Run tests
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
