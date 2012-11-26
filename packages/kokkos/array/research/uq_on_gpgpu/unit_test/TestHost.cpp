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

#include <iostream>

#include <KokkosArray_Host.hpp>

#include <KokkosArray_ProductTensor.hpp>
#include <KokkosArray_LegendrePolynomial.hpp>
#include <KokkosArray_SymmetricDiagonalSpec.hpp>
#include <KokkosArray_StochasticProductTensor.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>
#include <KokkosArray_CrsMatrix.hpp>


//

#include <Host/KokkosArray_Host_ProductTensor.hpp>
#include <Host/KokkosArray_Host_StochasticProductTensor.hpp>
#include <Host/KokkosArray_Host_SymmetricDiagonalSpec.hpp>
#include <Host/KokkosArray_Host_BlockCrsMatrix.hpp>
#include <Host/KokkosArray_Host_CrsMatrix.hpp>

//

#include <TestBlockCrsMatrix.hpp>
#include <TestTensorCrsMatrix.hpp>
#include <TestStochastic.hpp>

int mainHost(bool test_flat, bool test_orig, bool test_block)
{
  const size_t node_count = KokkosArray::Host::detect_node_count();
  const size_t node_thread_count = KokkosArray::Host::detect_node_core_count() / 2 ;

  KokkosArray::Host::initialize( node_count , node_thread_count );

//  unit_test::test_dense<KokkosArray::Host>();
//  unit_test::test_diagonal<KokkosArray::Host>();
//  unit_test::test_other<KokkosArray::Host>();

  unit_test::test_integration<1>();
  unit_test::test_integration<2>();
  unit_test::test_integration<3>();
  unit_test::test_integration<4>();
  unit_test::test_integration<5>();
  unit_test::test_integration<6>();
  unit_test::test_integration<7>();
  unit_test::test_integration<8>();
  unit_test::test_integration<9>();
  unit_test::test_integration<10>();
  unit_test::test_inner_product_legengre_polynomial<10,KokkosArray::Host>();
  unit_test::test_triple_product_legendre_polynomial<4,KokkosArray::Host>();

  unit_test::test_product_tensor<KokkosArray::Host,KokkosArray::SparseProductTensor>( std::vector<int>( 2 , 1 ) );
  unit_test::test_product_tensor<KokkosArray::Host,KokkosArray::SparseProductTensor>( std::vector<int>( 3 , 2 ) );
  unit_test::test_product_tensor<KokkosArray::Host,KokkosArray::SparseProductTensor>( std::vector<int>( 5 , 1 ) );

  unit_test::test_block_crs_matrix<KokkosArray::Host>( 1 , 2 );
  unit_test::test_block_crs_matrix<KokkosArray::Host>( 1 , 5 );
  unit_test::test_block_crs_matrix<KokkosArray::Host>( 2 , 1 );
  unit_test::test_block_crs_matrix<KokkosArray::Host>( 3 , 1 );

  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Host,long>( 1 , 2 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Host,long>( 1 , 5 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Host,long>( 2 , 1 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Host,long>( 5 , 1 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Host,long>( 5 , 5 );

  std::cout << "Stress tests:" << std::endl ;

  unit_test::test_block_crs_matrix<KokkosArray::Host>( 10 , 8 );
  unit_test::test_block_crs_matrix<KokkosArray::Host>( 11 , 8 );
  unit_test::test_block_crs_matrix<KokkosArray::Host>( 12 , 10 );
  unit_test::test_block_crs_matrix<KokkosArray::Host>( 13 , 10 );
  unit_test_tensor::test_tensor_crs_matrix<KokkosArray::Host,long>( 100 , 10 );

  std::cout << std::endl << "\"Host Performance with "
            << node_count * node_thread_count << " threads\"" << std::endl ;
  unit_test::performance_test_driver<KokkosArray::Host>(
    test_flat, test_orig, test_block);

  KokkosArray::Host::finalize();

  return 0 ;
}



