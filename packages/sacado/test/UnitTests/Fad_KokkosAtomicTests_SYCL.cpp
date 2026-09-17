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

// DFad is not tested on SYCL.  Every atomic test builds a local Fad value from
// the View -- "local_scalar_type x = m_v(i)" in AtomicKernel -- which for DFad
// is a device-side allocation, and SYCL has none:  neither operator new nor
// malloc resolves in device code.  SFad and SLFad keep their storage on the
// stack and are unaffected.
#define SACADO_TEST_DFAD 0

#include "Fad_KokkosAtomicTests.hpp"

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
