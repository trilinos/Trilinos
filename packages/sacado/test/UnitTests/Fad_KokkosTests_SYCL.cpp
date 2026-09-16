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

// Temporarily disable DFad testing on SYCL, matching the HIP backend.  DFad's
// temporary allocations need device-side "new", which SYCL does not provide.
#ifdef KOKKOS_ENABLE_SYCL
#define SACADO_TEST_DFAD 0
#else
#define SACADO_TEST_DFAD 1
#endif

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
