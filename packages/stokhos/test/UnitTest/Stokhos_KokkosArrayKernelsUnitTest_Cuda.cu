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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_KokkosKernelsUnitTest.hpp"

#include <Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_Cuda_ProductTensor.hpp>
#include <Cuda/Kokkos_Cuda_CrsProductTensorLegendre.hpp>
#include <Cuda/Kokkos_Cuda_StochasticProductTensor.hpp>
#include <Cuda/Kokkos_Cuda_SymmetricDiagonalSpec.hpp>
#include <Cuda/Kokkos_Cuda_CrsMatrix.hpp>
#include <Cuda/Kokkos_Cuda_BlockCrsMatrix.hpp>

using namespace KokkosKernelsUnitTest;

extern UnitTestSetup setup;

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsMatrixFree_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  bool test_block = true;
  
  success = test_crs_matrix_free<Scalar,Device>(setup, test_block, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsProductLegendre_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  
  success = test_crs_product_legendre<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsDenseBlock_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  
  success = test_crs_dense_block<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsFlatCommuted_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  
  success = test_crs_flat_commuted<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsFlatOriginal_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  
  success = test_crs_flat_original<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsProductTensor_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  
  success = test_crs_product_tensor<Scalar,Device,Kokkos::CrsProductTensor>(setup, out);
}
