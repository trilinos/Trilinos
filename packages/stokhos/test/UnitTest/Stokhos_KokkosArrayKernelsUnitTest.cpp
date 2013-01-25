// $Id$ 
// $Source$ 
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

#include "Stokhos_KokkosArrayKernelsUnitTest.hpp"

#include <Host/KokkosArray_Host_ProductTensor.hpp>
#include <Host/KokkosArray_Host_StochasticProductTensor.hpp>
#include <Host/KokkosArray_Host_SymmetricDiagonalSpec.hpp>
#include <Host/KokkosArray_Host_CrsMatrix.hpp>
#include <Host/KokkosArray_Host_BlockCrsMatrix.hpp>

#include <KokkosArray_Cuda.hpp>

using namespace KokkosArrayKernelsUnitTest;

UnitTestSetup setup;

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsMatrixFree_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  bool test_block = true;
  
  success = test_crs_matrix_free<Scalar,Device>(setup, test_block, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsProductLegendre_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_product_legendre<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsDenseBlock_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_dense_block<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsFlatCommuted_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_flat_commuted<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsFlatOriginal_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_flat_original<Scalar,Device>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, CrsProductTensor_Host ) {
  typedef double Scalar;
  typedef KokkosArray::Host Device;
  
  success = test_crs_product_tensor<Scalar,Device,KokkosArray::CrsProductTensor>(setup, out);
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, ProductLegendreCijk ) {
  success = true; 
  
  const std::vector<unsigned> var_degree( setup.d , setup.p );
  const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
    tensor( var_degree );
  
  typename UnitTestSetup::Cijk_type::k_iterator k_begin = 
    setup.Cijk->k_begin();
  typename UnitTestSetup::Cijk_type::k_iterator k_end = 
    setup.Cijk->k_end();
  for (typename UnitTestSetup::Cijk_type::k_iterator k_it=k_begin; 
       k_it!=k_end; ++k_it) {
    int k = index(k_it);
    typename UnitTestSetup::Cijk_type::kj_iterator j_begin = 
      setup.Cijk->j_begin(k_it);
    typename UnitTestSetup::Cijk_type::kj_iterator j_end = 
      setup.Cijk->j_end(k_it);
    for (typename UnitTestSetup::Cijk_type::kj_iterator j_it = j_begin; 
	 j_it != j_end; ++j_it) {
      int j = index(j_it);
      typename UnitTestSetup::Cijk_type::kji_iterator i_begin = 
	setup.Cijk->i_begin(j_it);
      typename UnitTestSetup::Cijk_type::kji_iterator i_end =  
	setup.Cijk->i_end(j_it);
      for (typename UnitTestSetup::Cijk_type::kji_iterator i_it = i_begin; 
	   i_it != i_end; 
	   ++i_it) {
	int i = index(i_it);
	double c = value(i_it);
	
	int ii = setup.perm[i];
	int jj = setup.perm[j];
	int kk = setup.perm[k];
	
	if (!tensor.is_non_zero(ii,jj,kk)) {
	  success = false;
	  out << "Cijk entry (" << ii << "," << jj << "," << kk 
	      << ") is zero!" << std::endl;
	}
	else {
	  double c2 = tensor(ii,jj,kk);
	  if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
	    out << "(" << i << "," << j << "," << k << "):  " << c 
		<< " == " << c2 << " failed!" << std::endl;
	    success = false;
	  }
	}
      }
    }
  }
}

TEUCHOS_UNIT_TEST( Stokhos_KokkosArrayKernels, ProductTensorCijk ) {
  success = true; 

  typedef double value_type;
  typedef KokkosArray::Host Device;
  typedef KokkosArray::NormalizedLegendrePolynomialBases<8> polynomial ;
  typedef KokkosArray::StochasticProductTensor< value_type , polynomial , Device , KokkosArray::CrsProductTensor > tensor_type ;
  
  const std::vector<unsigned> var_degree( setup.d , setup.p );

  tensor_type tensor = 
    KokkosArray::create_product_tensor< tensor_type >( var_degree );

  for (int i=0; i<setup.stoch_length; ++i) {
    const size_t iEntryBeg = tensor.tensor().entry_begin(i);
    const size_t iEntryEnd = tensor.tensor().entry_end(i);
    for (size_t iEntry = iEntryBeg ; iEntry < iEntryEnd ; ++iEntry ) {
      const size_t j = tensor.tensor().coord(iEntry,0);
      const size_t k = tensor.tensor().coord(iEntry,1);
      double c2 = tensor.tensor().value(iEntry);
      if (j == k) c2 *= 2.0;

      int ii = setup.inv_perm2[i];
      int jj = setup.inv_perm2[j];
      int kk = setup.inv_perm2[k];
      double c = setup.Cijk->getValue(ii,jj,kk);
      
      if (std::abs(c-c2) > std::abs(c)*setup.rel_tol + setup.abs_tol) {
	out << "(" << ii << "," << jj << "," << kk << "):  " << c 
	    << " == " << c2 << " failed!" << std::endl;
	success = false;
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
  KokkosArray::Host::initialize( gang_count , gang_worker_count );

  // Initialize Cuda
  KokkosArray::Cuda::initialize( KokkosArray::Cuda::SelectDevice(0) );

  // Run tests
  int ret = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  // Finish up
  KokkosArray::Host::finalize();
  KokkosArray::Cuda::finalize();

  return ret;
}
