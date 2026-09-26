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

// Re-test SYCL with hierarchical parallelism turned on (experimental), for
// SFad only.  With a fad size of 128 over a vector width of 32 each work item
// holds 4 derivative components, exercising the statically sized partitioned
// path.  Nothing here allocates.
#define SACADO_GPU_HIERARCHICAL 1

#include "Kokkos_Macros.hpp"

#define GLOBAL_FAD_SIZE 128

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

// Instantiate tests for SYCL device
using Kokkos::SYCL;
VIEW_FAD_TESTS_FDC(  SFadType , SYCL )
VIEW_FAD_TESTS_SFDC( SFadType , SYCL )

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
