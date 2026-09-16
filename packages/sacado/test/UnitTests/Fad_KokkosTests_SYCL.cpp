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

#include "Kokkos_Macros.hpp"

// DFad is tested on SYCL as it is on Cuda.  Inside a Kokkos::View the
// derivative array belongs to the View's allocation, so no device-side
// allocation is involved; the Fad temporaries a few kernels create fall back
// to operator new in device code, which the oneAPI/DPC++ extensions support.
#define SACADO_TEST_DFAD 1

#include "Fad_KokkosTests.hpp"

// Instantiate tests for SYCL device
using Kokkos::SYCL;
VIEW_FAD_TESTS_D( SYCL )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize SYCL
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize SYCL
  Kokkos::finalize();

  return res;
}
