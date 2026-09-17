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

// DFad on SYCL, in its own driver because it cannot run the whole test set.
//
// A DFad held in a Kokkos::View costs nothing extra on device:  the derivative
// array belongs to the View's allocation, made on the host, and the View hands
// out ViewFads that point into it.  That is what this driver covers.
//
// What SYCL cannot do is construct a DFad *value* in device code, which needs
// a device-side allocation.  Neither operator new nor malloc resolves in SYCL
// device code -- the kernel fails to build with "Unresolved Symbol" -- so the
// handful of tests that materialize a Fad value are compiled out here.  Unlike
// SFad and SLFad, whose storage is on the stack, DFad has no way around it.
#define SACADO_TEST_DEVICE_ALLOC 0

#define SACADO_TEST_DFAD 1

#include "Fad_KokkosTests.hpp"

// Instantiate tests for SYCL device, DFad only.  SFad and SLFad are covered by
// Fad_KokkosTests_SYCL, which runs the full set including the device-allocating
// tests.
using Kokkos::SYCL;
VIEW_FAD_TESTS_FD( DFadType, SYCL )

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
