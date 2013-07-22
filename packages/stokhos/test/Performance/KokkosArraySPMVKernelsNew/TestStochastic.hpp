// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
#include <utility>
#include <cmath>
#include <iostream>

#include "Kokkos_Host.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "Stokhos_CrsMatrix.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"
#include "Stokhos_CrsProductTensor.hpp"
#include "Stokhos_TiledCrsProductTensor.hpp"
#include "Stokhos_FlatSparse3Tensor.hpp"
#include "Stokhos_FlatSparse3Tensor_kji.hpp"
#include "Stokhos_LexicographicBlockSparse3Tensor.hpp"
#include "Stokhos_LinearSparse3Tensor.hpp"

#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_JacobiBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_TotalOrderBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_LTBSparse3Tensor.hpp"

namespace unit_test {

template< typename IntType >
inline
IntType map_fem_graph_coord( const IntType & N ,
                             const IntType & i ,
                             const IntType & j ,
                             const IntType & k )
{
  return k + N * ( j + N * i );
}

inline
size_t generate_fem_graph( size_t N ,
                           std::vector< std::vector<size_t> > & graph )
{
  graph.resize( N * N * N , std::vector<size_t>() );

  size_t total = 0 ;

  for ( int i = 0 ; i < (int) N ; ++i ) {
  for ( int j = 0 ; j < (int) N ; ++j ) {
  for ( int k = 0 ; k < (int) N ; ++k ) {

    const size_t row = map_fem_graph_coord((int)N,i,j,k);

    graph[row].reserve(27);

    for ( int ii = -1 ; ii < 2 ; ++ii ) {
    for ( int jj = -1 ; jj < 2 ; ++jj ) {
    for ( int kk = -1 ; kk < 2 ; ++kk ) {
      if ( 0 <= i + ii && i + ii < (int) N &&
           0 <= j + jj && j + jj < (int) N &&
           0 <= k + kk && k + kk < (int) N ) {
        size_t col = map_fem_graph_coord((int)N,i+ii,j+jj,k+kk);

        graph[row].push_back(col);
      }
    }}}
    total += graph[row].size();
  }}}

  return total ;
}

template <typename Scalar> struct ScalarTolerances {};

template <> struct ScalarTolerances<float> {
  typedef float scalar_type;
  static scalar_type sparse_cijk_tol() { return 1e-5; }
};

template <> struct ScalarTolerances<double> {
  typedef double scalar_type;
  static scalar_type sparse_cijk_tol() { return 1e-12; }
};

template< typename ScalarType , typename TensorType, class Device >
std::vector<double>
test_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type** ,
                             Kokkos::LayoutLeft ,
                             Device > block_vector_type ;

   typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

    typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::JacobiBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,order_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL);
  for (size_t i=0; i<num_KL; i++) {
    if (symmetric)
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,1.0,true));
    else
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,2.0,true));
  }
  RCP<const product_basis_type> basis =
    rcp(new product_basis_type(
          bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *basis,
                                                             *Cijk );
  matrix.graph = Kokkos::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();

  matrix.values = block_vector_type( "matrix" , inner_length , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM =
    Kokkos::create_mirror( matrix.values );

  for ( size_t i=0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1.0;
    }
  }

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1.0 ;
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.tensor().num_flops();
  const double flops = 1.0e-9*graph_length*flops_per_block;

  // std::cout << "tensor: flops = " << flops
  //        << " time = " << seconds_per_iter << std::endl;

  std::vector<double> perf(6) ;

  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops / seconds_per_iter;
  perf[3] = matrix.block.tensor().entry_count();
  perf[4] = inner_length ;
  perf[5] = flops_per_block;

  return perf ;
}

template< typename ScalarType , class Device >
std::vector<double>
test_tiled_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type** ,
                             Kokkos::LayoutLeft ,
                             Device > block_vector_type ;

  typedef Stokhos::TiledCrsProductTensor<ScalarType,Device> TensorType;
  typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

  typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::JacobiBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,order_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL);
  for (size_t i=0; i<num_KL; i++) {
    if (symmetric)
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,1.0,true));
    else
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,2.0,true));
  }
  RCP<const product_basis_type> basis =
    rcp(new product_basis_type(
          bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  Teuchos::ParameterList params;
  params.set("Tile Size", 128);
  params.set("Max Tiles", 10000);
  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *basis, *Cijk,
                                                             params);
  matrix.graph = Kokkos::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();

  matrix.values = block_vector_type( "matrix" , inner_length , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM =
    Kokkos::create_mirror( matrix.values );

  for ( size_t i=0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1.0;
    }
  }

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1.0 ;
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.tensor().num_flops();
  const double flops = 1.0e-9*graph_length*flops_per_block;

  // std::cout << "tensor: flops = " << flops
  //        << " time = " << seconds_per_iter << std::endl;

  std::vector<double> perf(6) ;

  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops / seconds_per_iter;
  perf[3] = matrix.block.tensor().entry_count();
  perf[4] = inner_length ;
  perf[5] = flops_per_block;

  return perf ;
}

template< typename ScalarType , class Device >
std::vector<double>
test_lexo_block_tensor(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type** ,
                             Kokkos::LayoutLeft ,
                             Device > block_vector_type ;

  typedef Stokhos::LexicographicBlockSparse3Tensor<value_type,Device> TensorType;
   typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

   typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::JacobiBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,order_type> product_basis_type;
  typedef Stokhos::LTBSparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL);
  for (size_t i=0; i<num_KL; i++) {
    if (symmetric)
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,1.0,true));
    else
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,2.0,true));
  }
  RCP<const product_basis_type> basis =
    rcp(new product_basis_type(
          bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  RCP<Cijk_type> Cijk =
    Stokhos::computeTripleProductTensorLTBBlockLeaf(*basis, symmetric);

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *basis,
                                                             *Cijk );
  matrix.graph = Kokkos::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();

  matrix.values = block_vector_type( "matrix" , inner_length , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM =
    Kokkos::create_mirror( matrix.values );

  for ( size_t i=0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1.0;
    }
  }

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1.0 ;
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.tensor().num_flops();
  const double flops = 1.0e-9*graph_length*flops_per_block;

  // std::cout << "tensor: flops = " << flops
  //        << " time = " << seconds_per_iter << std::endl;

  std::vector<double> perf(6) ;

  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops / seconds_per_iter;
  perf[3] = matrix.block.tensor().num_non_zeros();
  perf[4] = inner_length ;
  perf[5] = flops_per_block;

  return perf ;
}

template< typename ScalarType , class Device >
std::vector<double>
test_linear_tensor(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type** ,
                             Kokkos::LayoutLeft ,
                             Device > block_vector_type ;

  typedef Stokhos::LinearSparse3Tensor<value_type,Device,4> TensorType;
  typedef Stokhos::StochasticProductTensor< value_type , TensorType , Device > tensor_type ;

  typedef Stokhos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::JacobiBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,order_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL);
  for (size_t i=0; i<num_KL; i++) {
    if (symmetric)
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,1.0,true));
    else
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,2.0,true));
  }
  RCP<const product_basis_type> basis =
    rcp(new product_basis_type(
          bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  Teuchos::ParameterList params;
  params.set("Symmetric", symmetric);
  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *basis,
                                                             *Cijk,
                                                             params );
  matrix.graph = Kokkos::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length         = matrix.block.tensor().dimension();
  const size_t inner_length_aligned = matrix.block.tensor().aligned_dimension();

  matrix.values = block_vector_type( "matrix" , inner_length_aligned , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length_aligned , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length_aligned , outer_length );

  typename block_vector_type::HostMirror hM =
    Kokkos::create_mirror( matrix.values );

  for ( size_t i=0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1.0;
    }
  }

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1.0 ;
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.tensor().num_flops();
  const double flops = 1.0e-9*graph_length*flops_per_block;

  // std::cout << "tensor: flops = " << flops
  //        << " time = " << seconds_per_iter << std::endl;

  std::vector<double> perf(6) ;

  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops / seconds_per_iter;
  perf[3] = matrix.block.tensor().num_non_zeros();
  perf[4] = inner_length ;
  perf[5] = flops_per_block;

  return perf ;
}

template< typename ScalarType , class Device , class SparseMatOps >
std::vector<double>
test_original_matrix_free_vec(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool test_block ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::JacobiBasis<int,value_type> basis_type;
  typedef Stokhos::LexographicLess<Stokhos::MultiIndex<int> > order_type;
  typedef Stokhos::TotalOrderBasis<int,value_type,order_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL);
  for (size_t i=0; i<num_KL; i++) {
    if (symmetric)
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,1.0,true));
    else
      bases[i] = Teuchos::rcp(new basis_type(var_degree[i],1.0,2.0,true));
  }
  RCP<const product_basis_type> basis =
    rcp(new product_basis_type(
          bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  const size_t outer_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor();

  //------------------------------

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type ;
  typedef Kokkos::CrsArray<int,Device,void,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  typedef Kokkos::View<value_type[],Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  std::vector<vec_type> x( outer_length ) ;
  std::vector<vec_type> y( outer_length ) ;
  std::vector<vec_type> tmp( outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = Kokkos::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

    matrix[block].values = vec_type( "matrix" , graph_length );

    x[block]   = vec_type( "x" , inner_length );
    y[block]   = vec_type( "y" , inner_length );
    tmp[block] = vec_type( "tmp" , inner_length );

    typename vec_type::HostMirror hM =
      Kokkos::create_mirror( matrix[block].values );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1.0 ;
    }

    Kokkos::deep_copy( matrix[block].values , hM );

    typename vec_type::HostMirror hx =
      Kokkos::create_mirror( x[block] );
    typename vec_type::HostMirror hy =
      Kokkos::create_mirror( y[block] );

    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      hx(i) = 1.0 ;
      hy(i) = 0.0 ;
    }

    Kokkos::deep_copy( x[block] , hx );
    Kokkos::deep_copy( y[block] , hy );
  }

  SparseMatOps smo;
  Kokkos::Impl::Timer clock ;
  int n_apply = 0;
  int n_add = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    n_apply = 0;
    n_add = 0;
    typename Cijk_type::k_iterator k_begin = Cijk->k_begin();
    typename Cijk_type::k_iterator k_end = Cijk->k_end();
    for (typename Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int nj = Cijk->num_j(k_it);
      if (nj > 0) {
        int k = index(k_it);
        typename Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
        typename Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
        std::vector<vec_type> xx(nj), yy(nj);
        int jdx = 0;
        for (typename Cijk_type::kj_iterator j_it = j_begin; j_it != j_end;
             ++j_it) {
          int j = index(j_it);
          xx[jdx] = x[j];
          yy[jdx] = tmp[j];
          jdx++;
        }
        Stokhos::multiply( matrix[k] , xx , yy, test_block ,smo );
        n_apply += nj;
        jdx = 0;
        for (typename Cijk_type::kj_iterator j_it = j_begin; j_it != j_end;
             ++j_it) {
          typename Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
          typename Cijk_type::kji_iterator i_end =  Cijk->i_end(j_it);
          for (typename Cijk_type::kji_iterator i_it = i_begin; i_it != i_end;
               ++i_it) {
            int i = index(i_it);
            value_type c = value(i_it);
            Stokhos::update( value_type(1.0) , y[i] , c , yy[jdx] );
            ++n_add;
          }
          jdx++;
        }
      }
    }

  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 1.0e-9*(2.0*static_cast<double>(n_apply)*graph_length+
                               static_cast<double>(n_add)*inner_length);

  // std::cout << "mat-free: flops = " << flops
  //        << " time = " << seconds_per_iter << std::endl;

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length;
  perf[1] = seconds_per_iter ;
  perf[2] = flops/seconds_per_iter;
  perf[3] = flops;

  return perf;
}

template< class Scalar, class Device , class SparseMatOps >
void performance_test_driver_poly( const int pdeg ,
                                   const int minvar ,
                                   const int maxvar ,
                                   const int nGrid ,
                                   const int nIter ,
                                   const bool test_block ,
                                   const bool symmetric )
{
  // bool do_flat_sparse =
  //   Kokkos::Impl::is_same<Device,Kokkos::Host>::value ;

  std::cout.precision(8);

  //------------------------------

  std::vector< std::vector<size_t> > fem_graph ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );
  std::cout << std::endl << "\"FEM NNZ = " << graph_length << "\"" << std::endl;

  std::cout << std::endl
            << "\"#nGrid\" , "
            << "\"#Variable\" , "
            << "\"PolyDegree\" , "
            << "\"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
            << "\"Original-Matrix-Free-Block-MXV-Time\" , "
            << "\"Original-Matrix-Free-Block-MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block-MXV-GFLOPS\" , "
            << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
            // << "\"Block-Tiled Crs-Tensor MXV-Speedup\" , "
            // << "\"Block-Tiled Crs-3-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(
        var_degree , nGrid , nIter , symmetric );

    // const std::vector<double> perf_tiled_crs_tensor =
    //   test_tiled_product_tensor_matrix<Scalar,Device>(
    //     var_degree , nGrid , nIter , symmetric );

    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<Scalar,Device,SparseMatOps>(
        var_degree , nGrid , nIter , test_block , symmetric );

    std::cout << nGrid << " , "
              << nvar << " , "
              << pdeg << " , "
              << perf_crs_tensor[4] << " , "
              << perf_crs_tensor[3] << " , "
              << perf_original_mat_free_block[0] << " , "
              << perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
              << perf_original_mat_free_block[1] / perf_crs_tensor[1] << " , "
              // << perf_crs_tensor[1] << " , "
              << perf_crs_tensor[2] << " , "
              // << perf_original_mat_free_block[1] / perf_tiled_crs_tensor[1] << " , "
              // << perf_tiled_crs_tensor[2]
              << std::endl ;
  }

  //------------------------------
}

template< class Scalar, class Device , class SparseMatOps >
void performance_test_driver_poly_deg( const int nvar ,
                                       const int minp ,
                                       const int maxp ,
                                       const int nGrid ,
                                       const int nIter ,
                                       const bool test_block ,
                                       const bool symmetric )
{
  bool do_flat_sparse =
    Kokkos::Impl::is_same<Device,Kokkos::Host>::value ;

  std::cout.precision(8);

  //------------------------------

  std::vector< std::vector<size_t> > fem_graph ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );
  std::cout << std::endl << "\"FEM NNZ = " << graph_length << "\"" << std::endl;

  std::cout << std::endl
            << "\"#nGrid\" , "
            << "\"#Variable\" , "
            << "\"PolyDegree\" , "
            << "\"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
            << "\"Original-Matrix-Free-Block-MXV-Time\" , "
            << "\"Original-Matrix-Free-Block-MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block-MXV-GFLOPS\" , "
            << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , ";
  if (do_flat_sparse)
    std::cout << "\"Block-Lexicographic-Sparse-3-Tensor MXV-Speedup\" , "
              << "\"Block-Lexicographic-Sparse-3-Tensor MXV-GFLOPS\" , "
              << "\"Lexicographic FLOPS / Crs FLOPS\" , ";
  std::cout << std::endl ;

  for ( int p = minp ; p <= maxp ; ++p ) {
    std::vector<int> var_degree( nvar , p );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(
        var_degree , nGrid , nIter , symmetric );

    std::vector<double> perf_lexo_sparse_3_tensor;
    if (do_flat_sparse) {
      perf_lexo_sparse_3_tensor =
        test_lexo_block_tensor<Scalar,Device>( var_degree , nGrid , nIter , symmetric );
    }

    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<Scalar,Device,SparseMatOps>(
        var_degree , nGrid , nIter , test_block , symmetric );

    std::cout << nGrid << " , "
              << nvar << " , "
              << p << " , "
              << perf_crs_tensor[4] << " , "
              << perf_crs_tensor[3] << " , "
              << perf_original_mat_free_block[0] << " , "
              << perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
              << perf_original_mat_free_block[1] / perf_crs_tensor[1] << " , "
              << perf_crs_tensor[2] << " , ";
    if (do_flat_sparse) {
      std::cout << perf_original_mat_free_block[1] / perf_lexo_sparse_3_tensor[1] << " , "
                << perf_lexo_sparse_3_tensor[2] << " , "
                << perf_lexo_sparse_3_tensor[5] / perf_crs_tensor[5];
    }


    std::cout << std::endl ;
  }

  //------------------------------
}

template< class Scalar, class Device , class SparseMatOps >
void performance_test_driver_linear( const int minvar ,
                                     const int maxvar ,
                                     const int varinc ,
                                     const int nGrid ,
                                     const int nIter ,
                                     const bool test_block ,
                                     const bool symmetric )
{
  std::cout.precision(8);

  //------------------------------

  std::vector< std::vector<size_t> > fem_graph ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );
  std::cout << std::endl << "\"FEM NNZ = " << graph_length << "\"" << std::endl;

  std::cout << std::endl
            << "\"#nGrid\" , "
            << "\"#Variable\" , "
            << "\"PolyDegree\" , "
            << "\"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
            << "\"Original-Matrix-Free-Block-MXV-Time\" , "
            << "\"Original-Matrix-Free-Block-MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block-MXV-GFLOPS\" , "
            << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
            << "\"Linear-Sparse-3-Tensor MXV-Speedup\" , "
            << "\"Linear-Sparse-3-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; nvar+=varinc ) {
    std::vector<int> var_degree( nvar , 1 );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(
        var_degree , nGrid , nIter , symmetric );

    const std::vector<double> perf_linear_sparse_3_tensor =
      test_linear_tensor<Scalar,Device>( var_degree , nGrid , nIter , symmetric );

    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<Scalar,Device,SparseMatOps>(
        var_degree , nGrid , nIter , test_block , symmetric );

    std::cout << nGrid << " , "
              << nvar << " , "
              << 1 << " , "
              << perf_crs_tensor[4] << " , "
              << perf_crs_tensor[3] << " , "
              << perf_original_mat_free_block[0] << " , "
              << perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
              << perf_original_mat_free_block[1] / perf_crs_tensor[1] << " , "
              << perf_crs_tensor[2] << " , "
              << perf_original_mat_free_block[1] / perf_linear_sparse_3_tensor[1] << " , "
              << perf_linear_sparse_3_tensor[2] << " , "
              << std::endl ;
  }

  //------------------------------
}

template< class Scalar, class Device >
struct performance_test_driver {
  static void run(bool test_flat, bool test_orig, bool test_block, bool symmetric) {}
};

//----------------------------------------------------------------------------

}
