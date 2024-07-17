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

// Temporarily disable DFad testing on HIP. HIP does not support "new"
// on device so temporary allocations don't work.
#ifdef KOKKOS_ENABLE_HIP
#define SACADO_TEST_DFAD 0
#else
#define SACADO_TEST_DFAD 1
#endif

#include "Fad_KokkosAtomicTests.hpp"

// Instantiate tests for HIP device.
using Kokkos::HIP;
VIEW_FAD_TESTS_D( HIP )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize HIP
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize HIP
  Kokkos::finalize();

  return res;
}
