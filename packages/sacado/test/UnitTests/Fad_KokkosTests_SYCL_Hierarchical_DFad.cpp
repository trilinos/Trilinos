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

// Re-test SYCL with hierarchical DFad parallelism turned on (experimental).
//
// A DFad held in a Kokkos::View needs no device-side allocation:  the
// derivative array belongs to the View and the partitioned loop just strides
// over it.  That is what this covers.  The tests that build a Fad value in
// device code are compiled out, since those do allocate and SYCL has no
// device heap -- the memory pool the Cuda driver creates has no SYCL
// counterpart.
#define SACADO_GPU_HIERARCHICAL_DFAD 1

#include "Kokkos_Macros.hpp"

#define SACADO_TEST_DEVICE_ALLOC 0
#define SACADO_TEST_DFAD 1

#include "Fad_KokkosTests.hpp"

typedef Sacado::LayoutContiguous<Kokkos::LayoutLeft,32> LeftContiguous32;
typedef Sacado::LayoutContiguous<Kokkos::LayoutRight,32> RightContiguous32;
#undef VIEW_FAD_TESTS_FDC
#define VIEW_FAD_TESTS_FDC( F, D )                                      \
  VIEW_FAD_TESTS_FLD( F, LeftContiguous32, D )                          \
  VIEW_FAD_TESTS_FLD( F, RightContiguous32, D )

#undef VIEW_FAD_TESTS_SFDC
#define VIEW_FAD_TESTS_SFDC( F, D )                                     \
  VIEW_FAD_TESTS_SFLD( F, LeftContiguous32, D )                         \
  VIEW_FAD_TESTS_SFLD( F, RightContiguous32, D )

// Instantiate tests for SYCL device, DFad only
using Kokkos::SYCL;
VIEW_FAD_TESTS_FDC( DFadType , SYCL )

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
