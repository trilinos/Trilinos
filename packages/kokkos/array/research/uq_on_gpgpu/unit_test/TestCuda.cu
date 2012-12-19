/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <KokkosArray_Host.hpp>
#include <KokkosArray_Cuda.hpp>

#include <KokkosArray_ProductTensor.hpp>
#include <KokkosArray_LegendrePolynomial.hpp>
#include <KokkosArray_SymmetricDiagonalSpec.hpp>
#include <KokkosArray_StochasticProductTensor.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>
#include <KokkosArray_CrsMatrix.hpp>

//




#include <Host/KokkosArray_Host_ProductTensor.hpp>

#include <Cuda/KokkosArray_Cuda_SymmetricDiagonalSpec.hpp>
#include <Cuda/KokkosArray_Cuda_ProductTensor.hpp>
#include <Cuda/KokkosArray_Cuda_CrsProductTensorLegendre.hpp>
#include <Cuda/KokkosArray_Cuda_StochasticProductTensor.hpp>
#include <Cuda/KokkosArray_Cuda_BlockCrsMatrix.hpp>
#include <Cuda/KokkosArray_Cuda_CrsMatrix.hpp>

//

#include <TestBlockCrsMatrix.hpp>
#include <TestTensorCrsMatrix.hpp>
#include <TestStochastic.hpp>

namespace unit_test {

template<>
void performance_test_driver<KokkosArray::Cuda>(bool test_flat, bool test_orig, bool test_block)
{
  typedef KokkosArray::Cuda Device;

  int nGrid;
  int nIter; 
  bool print;

  // All methods compared against flat-original
  if (test_flat) {
    nGrid = 5 ;
    nIter = 1 ; 
    print = false ;
    performance_test_driver_all<Device>( 3 , 1 ,  9 , nGrid , nIter , print ,
    					 test_block );
    performance_test_driver_all<Device>( 5 , 1 ,  5 , nGrid , nIter , print ,
					 test_block );
  }

#ifdef HAVE_KOKKOSARRAY_STOKHOS
  // Just polynomial methods compared against original
  if (test_orig) {
    nGrid = 32 ;
    nIter = 1 ; 
    print = false ;
    performance_test_driver_poly<Device>( 3 , 1 , 12 , nGrid , nIter , print , 
					  test_block );
    performance_test_driver_poly<Device>( 5 , 1 ,  6 , nGrid , nIter , print ,
					  test_block );
  }
#endif

  //------------------------------

  /*
  std::cout << std::endl
            << "\"CRS flat-matrix ~27 nonzeros/row (CUDA uses cusparse)\""
            << std::endl
	    << "\"nGrid\" , "
            << "\"VectorSize\" , "
            << "\"MXV-Time\""
            << std::endl ;

  for ( int n_grid = 10 ; n_grid <= 100 ; n_grid += 5 ) {

    const std::pair<size_t,double> perf_flat =
      test_flat_matrix<double,Device>( n_grid , nIter , print );

    std::cout << n_grid << " , "
	      << perf_flat.first << " , "
              << perf_flat.second
              << std::endl ;
  }
  */

  //------------------------------
}

}

int mainCuda(bool test_flat, bool test_orig, bool test_block)
{
  typedef unsigned long long int IntType ;

  KokkosArray::Cuda::initialize();

//  unit_test::test_dense<KokkosArray::Cuda>();
//  unit_test::test_diagonal<KokkosArray::Cuda>();
//  unit_test::test_other<KokkosArray::Cuda>();

//  unit_test::test_inner_product_legengre_polynomial<10,KokkosArray::Cuda>();
//  unit_test::test_triple_product_legendre_polynomial<4,KokkosArray::Cuda>();

  unit_test::test_product_tensor<KokkosArray::Cuda,KokkosArray::SparseProductTensor>( std::vector<int>( 2 , 1 ) );
  unit_test::test_product_tensor<KokkosArray::Cuda,KokkosArray::SparseProductTensor>( std::vector<int>( 3 , 2 ) );
  unit_test::test_product_tensor<KokkosArray::Cuda,KokkosArray::SparseProductTensor>( std::vector<int>( 5 , 1 ) );

  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 1 , 2 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 1 , 5 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 2 , 1 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 2 , 2 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 3 , 1 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 3 , 2 );

  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Cuda,IntType>( 1 , 2 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Cuda,IntType>( 1 , 5 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Cuda,IntType>( 2 , 1 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Cuda,IntType>( 5 , 1 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Cuda,IntType>( 5 , 5 );

  std::cout << "Stress tests:" << std::endl ;

  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 10 , 8 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 11 , 8 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 12 , 10 );
  unit_test::test_block_crs_matrix<KokkosArray::Cuda>( 13 , 10 );

  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Cuda,IntType>( 100 , 10 );

  std::cout << std::endl << "\"Cuda Performance\"" << std::endl ;
  unit_test::performance_test_driver<KokkosArray::Cuda>(
    test_flat, test_orig, test_block);

  KokkosArray::Cuda::finalize();

  return 0 ;
}

