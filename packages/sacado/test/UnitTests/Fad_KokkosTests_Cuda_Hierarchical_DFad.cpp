// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Re-test cuda with hierarchical cuda parallelism turned on (experimental)
#define SACADO_VIEW_CUDA_HIERARCHICAL_DFAD 1
#define SACADO_KOKKOS_USE_MEMORY_POOL 1

#include "Fad_KokkosTests.hpp"
#include "Kokkos_Core.hpp"

typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft,32> LeftContiguous32;
typedef Kokkos::LayoutContiguous<Kokkos::LayoutRight,32> RightContiguous32;
#undef VIEW_FAD_TESTS_FDC
#define VIEW_FAD_TESTS_FDC( F, D )                                      \
  VIEW_FAD_TESTS_FLD( F, LeftContiguous32, D )                          \
  VIEW_FAD_TESTS_FLD( F, RightContiguous32, D )

#undef VIEW_FAD_TESTS_SFDC
#define VIEW_FAD_TESTS_SFDC( F, D )                                     \
  VIEW_FAD_TESTS_SFLD( F, LeftContiguous32, D )                         \
  VIEW_FAD_TESTS_SFLD( F, RightContiguous32, D )

// Instantiate tests for Cuda device
#if defined(KOKKOS_ENABLE_CUDA_UVM)
using Kokkos::Cuda;
VIEW_FAD_TESTS_FDC(  DFadType , Cuda )
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::InitArguments init_args;
  init_args.device_id = 0;
  Kokkos::initialize( init_args );
  Kokkos::print_configuration(std::cout);

#if defined(SACADO_KOKKOS_USE_MEMORY_POOL)
  Sacado::createGlobalMemoryPool(
    Kokkos::Cuda(),
    2*32*global_fad_size*global_num_rows*global_num_cols*sizeof(double),
    global_fad_size*sizeof(double),
    4*global_fad_size*sizeof(double),
    128*global_fad_size*sizeof(double));
#endif

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

#if defined(SACADO_KOKKOS_USE_MEMORY_POOL)
  Sacado::destroyGlobalMemoryPool(Kokkos::Cuda());
#endif

  // Finalize Cuda
  Kokkos::finalize();

  return res;
}
