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

#include "Stokhos_KokkosArrayKernelsUnitTestNew.hpp"

#include "Stokhos_Cuda_CrsMatrix.hpp"
#include "Stokhos_Cuda_BlockCrsMatrix.hpp"
#include "Stokhos_Cuda_StochasticProductTensor.hpp"
#include "Stokhos_Cuda_CrsProductTensor.hpp"
#include "Stokhos_Cuda_TiledCrsProductTensor.hpp"
#include "Stokhos_Cuda_LinearSparse3Tensor.hpp"

#include "Kokkos_Cuda.hpp"

using namespace KokkosKernelsUnitTest;

extern UnitTestSetup setup;

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsMatrixFree_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  typedef Stokhos::DefaultSparseMatOps SparseMatOps;
  bool test_block = true;

  success = test_crs_matrix_free<Scalar,Device,SparseMatOps>(
    setup, test_block, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, CrsProductTensor_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  typedef Stokhos::CrsProductTensor<Scalar,Device> Tensor;

  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, TiledCrsProductTensor_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  typedef Stokhos::TiledCrsProductTensor<Scalar,Device> Tensor;

  Teuchos::ParameterList params;
  params.set("Tile Size", 10);
  params.set("Max Tiles", 10000);
  success = test_crs_product_tensor<Scalar,Tensor,Device>(setup, out, params);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, LinearTensorSymmetric_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  const bool symmetric = true;

  UnitTestSetup s;
  s.setup(1, 40);
  success = test_linear_tensor<Scalar,Device,1>(s, out, symmetric);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosKernels, LinearTensorAsymmetric_Cuda ) {
  typedef double Scalar;
  typedef Kokkos::Cuda Device;
  const bool symmetric = false;

  UnitTestSetup s;
  s.setup(1, 40);
  success = test_linear_tensor<Scalar,Device,1>(s, out, symmetric);
}
