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

// Hierarchical SFad at a fad size large enough to put the kernels under real
// register pressure:  1024 components over a vector width of 32 is 32 doubles
// held per thread, against 4 for the 128-component driver.
//
// This is the correctness counterpart to the register-pressure sweep in
// SyclIndexProbe.  That sweep measures the vector width Kokkos delivers as a
// kernel grows; this checks whether Sacado still computes the right
// derivatives at that width.  Kokkos reduces the vector dimension to a
// kernel's own maximum sub-group size, which falls as registers get tight,
// while Sacado sizes the partitioned Fad type at compile time from the layout
// stride.  If the two disagree each thread writes past the end of its local
// Fad.  See SacadoSyclStrideIssue.txt.
#define SACADO_VIEW_CUDA_HIERARCHICAL 1

#include "Kokkos_Macros.hpp"

#define GLOBAL_FAD_SIZE 1024
#define SACADO_TEST_DFAD 0

#include "Fad_KokkosTests.hpp"

typedef Sacado::LayoutContiguous<Kokkos::LayoutLeft,32> LeftContiguous32;
typedef Sacado::LayoutContiguous<Kokkos::LayoutRight,32> RightContiguous32;

// Only the tests that launch a partitioned kernel over a small view.  The full
// set is deliberately not instantiated:  Rank8 allocates 100x1x2x3x4x5x6
// elements, which at this fad size would be roughly 590 MB on the device plus
// a host mirror.
#define LARGE_SFAD_TESTS( L, D )                                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ScalarAssign,   SFadType, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ValueAssign,    SFadType, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Multiply,       SFadType, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, MultiplyUpdate, SFadType, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, MultiplyConst,  SFadType, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicAdd,      SFadType, L, D )

using Kokkos::SYCL;
LARGE_SFAD_TESTS( LeftContiguous32,  SYCL )
LARGE_SFAD_TESTS( RightContiguous32, SYCL )

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
