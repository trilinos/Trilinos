// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

// Utilities
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Device
#include "Kokkos_Core.hpp"

// Kernels
#include "Stokhos_Cuda_CrsMatrix.hpp"
#include "Stokhos_Cuda_BlockCrsMatrix.hpp"
#include "Stokhos_Cuda_StochasticProductTensor.hpp"
#include "Stokhos_Cuda_SymmetricDiagonalSpec.hpp"
#include "Stokhos_Cuda_CrsProductTensor.hpp"
#include "Stokhos_Cuda_TiledCrsProductTensor.hpp"
#include "Stokhos_Cuda_SimpleTiledCrsProductTensor.hpp"
#include "Stokhos_Cuda_CooProductTensor.hpp"
#include "Stokhos_Cuda_LinearSparse3Tensor.hpp"

// Tests
#include "Stokhos_KokkosArrayKernelsUnitTest.hpp"

using namespace KokkosKernelsUnitTest;

UnitTestSetup<Kokkos::Cuda> setup;

// Test declarations
#include "Stokhos_KokkosArrayKernelsUnitTestDecl.hpp"

// Not all of the generic kernels work with Cuda because of the way
// the Cuda StochasticProductTensor specialization is done

using Kokkos::Cuda;

// The Tiled-Crs kernel seems to fail when OpenMP is the host device
#ifdef KOKKOS_ENABLE_OPENMP
#define TILED_CRS_TEST(SCALAR, DEVICE)
#else
#define TILED_CRS_TEST(SCALAR, DEVICE)                                         \
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, TiledCrsProductTensor, SCALAR, DEVICE )
#endif

#define UNIT_TEST_GROUP_SCALAR_CUDA( SCALAR ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFree, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFreeView, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFreeKokkos, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsMatrixFreeSingleCol, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsDenseBlock, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsFlatCommuted, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsFlatOriginal, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CrsProductTensor, SCALAR, Cuda ) \
  TILED_CRS_TEST(SCALAR, Cuda )                                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CooProductTensorPacked, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, CooProductTensorUnpacked, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, LinearTensorSymmetric, SCALAR, Cuda ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, LinearTensorAsymmetric, SCALAR, Cuda )

// Commenting this one out -- it may be generating memory errors
// TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Kokkos_SG_SpMv, SimpleTiledCrsProductTensor, SCALAR, Cuda )

UNIT_TEST_GROUP_SCALAR_CUDA(double)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // Initialize Cuda
  Kokkos::InitArguments init_args;
  init_args.device_id = 0;
  Kokkos::initialize( init_args );
  Kokkos::print_configuration( std::cout );

  // Setup (has to happen after initialization)
  setup.setup();

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  Kokkos::finalize();

  return ret;
}
