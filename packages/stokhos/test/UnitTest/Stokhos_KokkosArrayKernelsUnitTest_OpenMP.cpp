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
#include "Stokhos_ConfigDefs.h"
#include "Stokhos_OpenMP_CrsProductTensor.hpp"
#ifdef HAVE_STOKHOS_MKL
#include "Stokhos_OpenMP_MKL_CrsMatrix.hpp"
#endif

// Tests
#include "Stokhos_KokkosArrayKernelsUnitTest.hpp"

using namespace KokkosKernelsUnitTest;

UnitTestSetup<Kokkos::OpenMP> setup;

// Test declarations
#include "Stokhos_KokkosArrayKernelsUnitTestDecl.hpp"
#include "Stokhos_KokkosArrayKernelsUnitTest_Host.hpp"

// Tests using OpenMP device
using Kokkos::OpenMP;
UNIT_TEST_GROUP_SCALAR_DEVICE( double, OpenMP )
UNIT_TEST_GROUP_SCALAR_HOST_DEVICE( double, OpenMP )

#ifdef HAVE_STOKHOS_MKL
TEUCHOS_UNIT_TEST( Kokkos_SG_SpMv, double_OpenMP_CrsMatrixFree_MKL ) {
  typedef double Scalar;
  typedef Kokkos::OpenMP Device;
  typedef Stokhos::MKLMultiply SparseMatOps;
  success = test_crs_matrix_free<Scalar,Device,SparseMatOps>(
    setup, out);
}
#endif

int main( int argc, char* argv[] ) {
  // Setup the MPI session
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Kokkos
  Kokkos::initialize(argc, argv);

  // Setup (has to happen after initialization)
  setup.setup();

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
