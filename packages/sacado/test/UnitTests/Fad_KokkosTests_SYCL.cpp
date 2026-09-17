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

// DFad is covered by its own driver, Fad_KokkosTests_DFad_SYCL:  it cannot run
// the tests that construct a Fad value in device code, because SYCL provides no
// device-side allocation.  Keeping it out of here lets SFad and SLFad run the
// full set.
#define SACADO_TEST_DFAD 0

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
