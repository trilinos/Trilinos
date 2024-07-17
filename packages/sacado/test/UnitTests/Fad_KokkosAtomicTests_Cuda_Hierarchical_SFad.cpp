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

// Re-test cuda with hierarchical cuda parallelism turned on (experimental)
#define SACADO_VIEW_CUDA_HIERARCHICAL 1

#define GLOBAL_FAD_SIZE 64

#include "Fad_KokkosAtomicTests.hpp"

typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft,32> LeftContiguous32;
typedef Kokkos::LayoutContiguous<Kokkos::LayoutRight,32> RightContiguous32;
#undef VIEW_FAD_TESTS_FD
#define VIEW_FAD_TESTS_FD( F, D )                                       \
  VIEW_FAD_TESTS_FLD( F, LeftContiguous32, D )                          \
  VIEW_FAD_TESTS_FLD( F, RightContiguous32, D )

// Instantiate tests for Cuda device
#if SACADO_ENABLE_NEW_DESIGN
using Kokkos::Cuda;
VIEW_FAD_TESTS_FD( SFadType , Cuda )
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(0);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finalize Cuda
  Kokkos::finalize();

  return res;
}
