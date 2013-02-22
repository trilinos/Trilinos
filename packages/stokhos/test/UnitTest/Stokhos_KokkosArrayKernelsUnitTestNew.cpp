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

#include "Stokhos_Host_CrsMatrix.hpp"
#include "Stokhos_Host_BlockCrsMatrix.hpp"
#include "Stokhos_Host_StochasticProductTensor.hpp"
#include "Stokhos_Host_CrsProductTensor.hpp"
#include "Stokhos_Host_FlatSparse3Tensor.hpp"
#include "Stokhos_Host_FlatSparse3Tensor_kji.hpp"

#include "KokkosArray_Cuda.hpp"

using namespace KokkosArrayKernelsUnitTest;

UnitTestSetup setup;

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsMatrixFree_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  bool test_block = true;
  
  success = test_crs_matrix_free<Scalar,Device>(setup, test_block, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsProductTensor_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_product_tensor<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparse3Tensor_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_product_tensor<Scalar,Stokhos::FlatSparse3Tensor<Scalar,Device>,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparse3Tensor_kji_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_product_tensor<Scalar,Stokhos::FlatSparse3Tensor_kji<Scalar,Device>,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, ProductTensorCijk ) {
  success = true; 

  typedef double value_type;
  typedef KokkosArray::Host Device;
  typedef Stokhos::CrsProductTensor< value_type , Device > tensor_type ;

  tensor_type tensor = 
   Stokhos::create_product_tensor<Device>( *setup.basis, *setup.Cijk );

  for (int i=0; i<setup.stoch_length; ++i) {
    const size_t iEntryBeg = tensor.entry_begin(i);
    const size_t iEntryEnd = tensor.entry_end(i);
    for (size_t iEntry = iEntryBeg ; iEntry < iEntryEnd ; ++iEntry ) {
      const size_t j = tensor.coord(iEntry,0);
      const size_t k = tensor.coord(iEntry,1);
      double c2 = tensor.value(iEntry);
      if (j == k) c2 *= 2.0;

      int ii = setup.inv_perm[i];
      int jj = setup.inv_perm[j];
      int kk = setup.inv_perm[k];
      double c = setup.Cijk->getValue(ii,jj,kk);
      
      if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
	out << "(" << ii << "," << jj << "," << kk << "):  " << c 
	    << " == " << c2 << " failed!" << std::endl;
	success = false;
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparseCijk ) {
  success = true; 

  typedef double value_type;
  typedef KokkosArray::Host Device;
  typedef Stokhos::FlatSparse3Tensor< value_type , Device > tensor_type ;
  typedef size_t size_type;

  tensor_type tensor = 
   Stokhos::create_flat_sparse_3_tensor<Device>( *setup.basis, *setup.Cijk );

  for (int i=0; i<setup.stoch_length; ++i) {
    const size_type nk = tensor.num_k(i);
    const size_type kBeg = tensor.k_begin(i);
    const size_type kEnd = kBeg + nk;
    for (size_type kEntry = kBeg; kEntry < kEnd; ++kEntry) {
      const size_type k = tensor.k_coord(kEntry);
      const size_type nj = tensor.num_j(kEntry);
      const size_type jBeg = tensor.j_begin(kEntry);
      const size_type jEnd = jBeg + nj;
      for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
	const size_type j = tensor.j_coord(jEntry);
	double c2 = tensor.value(jEntry);
	if (j == k) c2 *= 2.0;
	double c = setup.Cijk->getValue(i,j,k);
	if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
	  out << "(" << i << "," << j << "," << k << "):  " << c 
	      << " == " << c2 << " failed!" << std::endl;
	  success = false;
	}
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, FlatSparseCijk_kji ) {
  success = true; 

  typedef double value_type;
  typedef KokkosArray::Host Device;
  typedef Stokhos::FlatSparse3Tensor_kji< value_type , Device > tensor_type ;
  typedef size_t size_type;

  tensor_type tensor = 
   Stokhos::create_flat_sparse_3_tensor_kji<Device>(*setup.basis, *setup.Cijk);
  const size_type nk = tensor.num_k();

  for ( size_type k = 0; k < nk; ++k) {
    const size_type nj = tensor.num_j(k);
    const size_type jBeg = tensor.j_begin(k);
    const size_type jEnd = jBeg + nj;
    for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
      const size_type j = tensor.j_coord(jEntry);
      const size_type ni = tensor.num_i(jEntry);
      const size_type iBeg = tensor.i_begin(jEntry);
      const size_type iEnd = iBeg + ni;
      for (size_type iEntry = iBeg; iEntry < iEnd; ++iEntry) {
	const size_type i = tensor.i_coord(iEntry);
	double c2 = tensor.value(iEntry);
	if (j == k) c2 *= 2.0;
	double c = setup.Cijk->getValue(i,j,k);
	if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
	  out << "(" << i << "," << j << "," << k << "):  " << c 
	      << " == " << c2 << " failed!" << std::endl;
	  success = false;
	}
      }
    }
  }
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  setup.setup();

  // Initialize host
  const size_t gang_count = 
    KokkosArray::Host::detect_gang_capacity();
  const size_t gang_worker_count = 
    KokkosArray::Host::detect_gang_worker_capacity() / 2 ;
  // const size_t gang_count = 1;
  // const size_t gang_worker_count = 1;
  KokkosArray::Host::initialize( gang_count , gang_worker_count );

#ifdef HAVE_KOKKOSARRAY_CUDA
  // Initialize Cuda
  KokkosArray::Cuda::initialize( KokkosArray::Cuda::SelectDevice(0) );
#endif

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  KokkosArray::Host::finalize();
#ifdef HAVE_KOKKOSARRAY_CUDA
  KokkosArray::Cuda::finalize();
#endif

  return ret;
}
