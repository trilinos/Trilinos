// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <utility>
#include <cmath>
#include <iostream>

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"

#include "Stokhos_Update.hpp"
#include "Stokhos_CrsMatrix.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"
#include "Stokhos_SymmetricDiagonalSpec.hpp"
#include "Stokhos_CrsProductTensor.hpp"
#include "Stokhos_TiledCrsProductTensor.hpp"
#include "Stokhos_SimpleTiledCrsProductTensor.hpp"
#include "Stokhos_CooProductTensor.hpp"
#include "Stokhos_FlatSparse3Tensor.hpp"
#include "Stokhos_FlatSparse3Tensor_kji.hpp"
#include "Stokhos_LexicographicBlockSparse3Tensor.hpp"
#include "Stokhos_LinearSparse3Tensor.hpp"

#if defined(HAVE_MPI) && 0
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_JacobiBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_TotalOrderBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Stokhos_LTBSparse3Tensor.hpp"
#include "Stokhos_Sparse3TensorUtilities.hpp"

#ifdef HAVE_STOKHOS_KOKKOSLINALG
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosBlas1_update.hpp"
#endif

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
  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length         = matrix.block.dimension();
  const size_t inner_length_aligned = matrix.block.aligned_dimension();

  matrix.values =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), inner_length_aligned , graph_length );

  block_vector_type x =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), inner_length_aligned , outer_length );
  block_vector_type y =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), inner_length_aligned , outer_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //------------------------------
  // Generate input multivector:

  Kokkos::deep_copy( x , ScalarType(1.0) );
  block_vector_type x0 =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("x"),
                       inner_length_aligned , outer_length );
  Kokkos::deep_copy( x0 , ScalarType(1.0) );

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Kokkos::deep_copy( x, x0 ); // akin to import
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.tensor().num_flops();
  const double flops = 1.0e-9*graph_length*flops_per_block;

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
test_product_tensor_diagonal_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type**,
                             Kokkos::LayoutLeft ,
                             Device > block_vector_type ;

  //------------------------------

  typedef Stokhos::BlockCrsMatrix< Stokhos::SymmetricDiagonalSpec< Device > ,
                                  value_type , Device > matrix_type ;

  typedef typename matrix_type::graph_type  graph_type ;

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
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;
  const size_t fem_graph_length = unit_test::generate_fem_graph( nGrid , fem_graph );

  const size_t stoch_length = basis->size();

  //------------------------------

  block_vector_type x = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), stoch_length , fem_length );
  block_vector_type y = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), stoch_length , fem_length );

  Kokkos::deep_copy( x , ScalarType(1.0) );

  //------------------------------
  // Generate CRS matrix of blocks with symmetric diagonal storage

  matrix_type matrix ;

  matrix.block  = Stokhos::SymmetricDiagonalSpec< Device >( stoch_length );
  matrix.graph  = Kokkos::create_staticcrsgraph<graph_type>(
    std::string("test product tensor graph") , fem_graph );
  matrix.values = block_vector_type(
    Kokkos::ViewAllocateWithoutInitializing("matrix"), matrix.block.matrix_size() , fem_graph_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = 2.0*stoch_length*stoch_length;
  const double flops = 1e-9*fem_graph_length*flops_per_block;

  std::vector<double> perf(6);
  perf[0] = fem_length * stoch_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops / seconds_per_iter;
  perf[3] = Cijk->num_entries();
  perf[4] = stoch_length;
  perf[5] = flops_per_block;
  return perf;
}

//----------------------------------------------------------------------------
// Flatten to a plain CRS matrix
//
//  Outer DOF == fem
//  Inner DOF == stochastic

template< typename ScalarType , class Device >
std::vector<double>
test_product_flat_commuted_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type* , Device > vector_type ;

  //------------------------------

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;

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
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;

  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Stochastic graph

  const size_t stoch_length = basis->size();
  std::vector< std::vector< int > > stoch_graph( stoch_length );
#if defined(HAVE_MPI) && 0
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
  Teuchos::RCP<Epetra_CrsGraph> cijk_graph = Stokhos::sparse3Tensor2CrsGraph(
    *basis, *Cijk, comm);
  for ( size_t i = 0 ; i < stoch_length ; ++i ) {
    int len = cijk_graph->NumGlobalIndices(i);
    stoch_graph[i].resize(len);
    int len2;
    cijk_graph->ExtractGlobalRowCopy(i, len, len2, &stoch_graph[i][0]);
  }

  //------------------------------
  // Generate flattened graph with FEM outer and stochastic inner

  const size_t flat_length = fem_length * stoch_length ;

  std::vector< std::vector<size_t> > flat_graph( flat_length );

  for ( size_t iOuterRow = 0 ; iOuterRow < fem_length ; ++iOuterRow ) {

    const size_t iOuterRowNZ = fem_graph[iOuterRow].size();

    for ( size_t iInnerRow = 0 ; iInnerRow < stoch_length ; ++iInnerRow ) {

      const size_t iInnerRowNZ = stoch_graph[ iInnerRow ].size(); ;
      const size_t iFlatRowNZ  = iOuterRowNZ * iInnerRowNZ ;
      const size_t iFlatRow    = iInnerRow + iOuterRow * stoch_length ;

      flat_graph[iFlatRow].resize( iFlatRowNZ );

      size_t iFlatEntry = 0 ;

      for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterRowNZ ; ++iOuterEntry ) {

        const size_t iOuterCol = fem_graph[iOuterRow][iOuterEntry];

        for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerRowNZ ; ++iInnerEntry ) {

          const size_t iInnerCol   = stoch_graph[iInnerRow][iInnerEntry] ;
          const size_t iFlatColumn = iInnerCol + iOuterCol * stoch_length ;

          flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;

          ++iFlatEntry ;
        }
      }
    }
  }

  //------------------------------

  vector_type x = vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), flat_length );
  vector_type y = vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), flat_length );

  Kokkos::deep_copy( x , ScalarType(1.0) );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = Kokkos::create_staticcrsgraph<matrix_graph_type>(
    std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.extent(0);

  matrix.values = matrix_values_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), flat_graph_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //Kokkos::write_matrix_market(matrix, "flat_commuted.mm");

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*flat_graph_length / seconds_per_iter;

  std::vector<double> perf(4);
  perf[0] = flat_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  perf[3] = flat_graph_length ;
  return perf;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Flatten to a plain CRS matrix
//
//  Outer DOF == stochastic
//  Inner DOF == fem

template< typename ScalarType , class Device >
std::vector<double>
test_product_flat_original_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type* , Device > vector_type ;

  //------------------------------

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;

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
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;

  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Stochastic graph

  const size_t stoch_length = basis->size();
  std::vector< std::vector< int > > stoch_graph( stoch_length );
#if defined(HAVE_MPI) && 0
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
  Teuchos::RCP<Epetra_CrsGraph> cijk_graph = Stokhos::sparse3Tensor2CrsGraph(
    *basis, *Cijk, comm);
  for ( size_t i = 0 ; i < stoch_length ; ++i ) {
    int len = cijk_graph->NumGlobalIndices(i);
    stoch_graph[i].resize(len);
    int len2;
    cijk_graph->ExtractGlobalRowCopy(i, len, len2, &stoch_graph[i][0]);
  }

  //------------------------------
  // Generate flattened graph with stochastic outer and FEM inner

  const size_t flat_length = fem_length * stoch_length ;

  std::vector< std::vector<size_t> > flat_graph( flat_length );

  for ( size_t iOuterRow = 0 ; iOuterRow < stoch_length ; ++iOuterRow ) {

    const size_t iOuterRowNZ = stoch_graph[iOuterRow].size();

    for ( size_t iInnerRow = 0 ; iInnerRow < fem_length ; ++iInnerRow ) {

      const size_t iInnerRowNZ = fem_graph[iInnerRow].size();
      const size_t iFlatRowNZ  = iOuterRowNZ * iInnerRowNZ ;
      const size_t iFlatRow    = iInnerRow + iOuterRow * fem_length ;

      flat_graph[iFlatRow].resize( iFlatRowNZ );

      size_t iFlatEntry = 0 ;

      for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterRowNZ ; ++iOuterEntry ) {

        const size_t iOuterCol = stoch_graph[ iOuterRow ][ iOuterEntry ];

        for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerRowNZ ; ++iInnerEntry ) {

          const size_t iInnerCol   = fem_graph[ iInnerRow][iInnerEntry];
          const size_t iFlatColumn = iInnerCol + iOuterCol * fem_length ;

          flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;
          ++iFlatEntry ;
        }
      }
    }
  }

  //------------------------------

  vector_type x = vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), flat_length );
  vector_type y = vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), flat_length );

  Kokkos::deep_copy( x , ScalarType(1.0) );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = Kokkos::create_staticcrsgraph<matrix_graph_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.extent(0);

  matrix.values = matrix_values_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), flat_graph_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //Kokkos::write_matrix_market(matrix, "flat_original.mm");

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*flat_graph_length / seconds_per_iter;

  std::vector<double> perf(4);
  perf[0] = flat_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  perf[3] = flat_graph_length ;
  return perf;
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
  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();
  const size_t inner_length_aligned = matrix.block.aligned_dimension();

  matrix.values =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), inner_length_aligned , graph_length );

  block_vector_type x =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), inner_length_aligned , outer_length );
  block_vector_type y =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), inner_length_aligned , outer_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //------------------------------
  // Generate input multivector:

  Kokkos::deep_copy( x , ScalarType(1.0) );

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

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
test_simple_tiled_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool symmetric )
{
  typedef ScalarType value_type ;
  typedef Kokkos::View< value_type** ,
                             Kokkos::LayoutLeft ,
                             Device > block_vector_type ;

  typedef Stokhos::SimpleTiledCrsProductTensor<ScalarType,Device> TensorType;
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
  matrix.block =
    Stokhos::create_stochastic_product_tensor< TensorType >( *basis, *Cijk,
                                                             params);
  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();
  const size_t inner_length_aligned = matrix.block.aligned_dimension();

  matrix.values =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), inner_length_aligned , graph_length );

  block_vector_type x =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), inner_length_aligned , outer_length );
  block_vector_type y =
    block_vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), inner_length_aligned , outer_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //------------------------------
  // Generate input multivector:

  Kokkos::deep_copy( x , ScalarType(1.0) );

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

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
  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();

  matrix.values = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), inner_length , graph_length );

  block_vector_type x = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), inner_length , outer_length );
  block_vector_type y = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), inner_length , outer_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //------------------------------
  // Generate input multivector:

  Kokkos::deep_copy( x , ScalarType(1.0) );

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

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
  matrix.graph = Kokkos::create_staticcrsgraph<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length         = matrix.block.tensor().dimension();
  const size_t inner_length_aligned = matrix.block.tensor().aligned_dimension();

  matrix.values = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), inner_length_aligned , graph_length );

  block_vector_type x = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("x"), inner_length_aligned , outer_length );
  block_vector_type y = block_vector_type( Kokkos::ViewAllocateWithoutInitializing("y"), inner_length_aligned , outer_length );

  Kokkos::deep_copy( matrix.values , ScalarType(1.0) );

  //------------------------------
  // Generate input multivector:

  Kokkos::deep_copy( x , ScalarType(1.0) );

  //------------------------------

  Device().fence();
  Kokkos::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Stokhos::multiply( matrix , x , y );
  }
  Device().fence();

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

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  typedef Kokkos::View<value_type*,Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  std::vector<vec_type> x( outer_length ) ;
  std::vector<vec_type> y( outer_length ) ;
  std::vector<vec_type> tmp( outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = Kokkos::create_staticcrsgraph<matrix_graph_type>( std::string("testing") , fem_graph );

    matrix[block].values = matrix_values_type( Kokkos::ViewAllocateWithoutInitializing("matrix"), graph_length );

    x[block]   = vec_type( Kokkos::ViewAllocateWithoutInitializing("x"), inner_length );
    y[block]   = vec_type( Kokkos::ViewAllocateWithoutInitializing("y"), inner_length );
    tmp[block] = vec_type( Kokkos::ViewAllocateWithoutInitializing("tmp"), inner_length );

    Kokkos::deep_copy( matrix[block].values , ScalarType(1.0) );

    Kokkos::deep_copy( x[block] , ScalarType(1.0) );
    Kokkos::deep_copy( y[block] , ScalarType(1.0) );
  }

  Device().fence();
  SparseMatOps smo;
  Kokkos::Timer clock ;
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
        Stokhos::multiply( matrix[k] , xx , yy, smo );
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
  Device().fence();

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

template< typename ScalarType , class Device , class SparseMatOps >
std::vector<double>
test_original_matrix_free_view(
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

  typedef Stokhos::CrsMatrix<value_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::graph_type matrix_graph_type;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  typedef Kokkos::View<value_type*, Kokkos::LayoutLeft, Device, Kokkos::MemoryUnmanaged> vec_type ;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Device> multi_vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  multi_vec_type x( Kokkos::ViewAllocateWithoutInitializing("x"),
                    inner_length, outer_length  ) ;
  multi_vec_type y("y", inner_length, outer_length ) ;
  multi_vec_type tmp_x( "tmp_x", inner_length, outer_length ) ;
  multi_vec_type tmp_y( "tmp_y", inner_length, outer_length ) ;

  Kokkos::deep_copy( x , ScalarType(1.0) );

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("testing") , fem_graph );

    matrix[block].values = matrix_values_type( "matrix" , graph_length );

    Kokkos::deep_copy( matrix[block].values , ScalarType(1.0) );
  }

  Device().fence();
  SparseMatOps smo;
  Kokkos::Timer clock ;
  int n_apply = 0;
  int n_add = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    typedef typename Cijk_type::k_iterator k_iterator;
    typedef typename Cijk_type::kj_iterator kj_iterator;
    typedef typename Cijk_type::kji_iterator kji_iterator;
    n_apply = 0;
    n_add = 0;
    k_iterator k_begin = Cijk->k_begin();
    k_iterator k_end = Cijk->k_end();
    for (k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      unsigned nj = Cijk->num_j(k_it);
      if (nj > 0) {
        int k = index(k_it);
        kj_iterator j_begin = Cijk->j_begin(k_it);
        kj_iterator j_end = Cijk->j_end(k_it);
        std::vector<int> j_indices(nj);
        unsigned jdx = 0;
        for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          int j = index(j_it);
          vec_type xx = Kokkos::subview( x, Kokkos::ALL(), j );
          vec_type tt = Kokkos::subview( tmp_x, Kokkos::ALL(), jdx++ );
          Kokkos::deep_copy(tt, xx);
        }
        multi_vec_type tmp_x_view =
          Kokkos::subview( tmp_x, Kokkos::ALL(),
                                           std::make_pair(0u,nj));
        multi_vec_type tmp_y_view =
          Kokkos::subview( tmp_y, Kokkos::ALL(),
                                           std::make_pair(0u,nj));
        Stokhos::multiply( matrix[k] , tmp_x_view , tmp_y_view, smo );
        n_apply += nj;
        jdx = 0;
        for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          vec_type tmp_y_view =
            Kokkos::subview( tmp_y, Kokkos::ALL(), jdx++ );
          kji_iterator i_begin = Cijk->i_begin(j_it);
          kji_iterator i_end =  Cijk->i_end(j_it);
          for (kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
            int i = index(i_it);
            value_type c = value(i_it);
            vec_type y_view = Kokkos::subview( y, Kokkos::ALL(), i );
            Stokhos::update( value_type(1.0) , y_view , c , tmp_y_view );
            ++n_add;
          }
        }
      }
    }

  }
  Device().fence();

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

#ifdef HAVE_STOKHOS_KOKKOSLINALG
template< typename ScalarType , class Device >
std::vector<double>
test_original_matrix_free_kokkos(
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

  typedef int ordinal_type;
  typedef KokkosSparse::CrsMatrix<value_type,ordinal_type,Device> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  typedef Kokkos::View<value_type*,Kokkos::LayoutLeft,Device, Kokkos::MemoryUnmanaged> vec_type ;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Device> multi_vec_type;

  std::vector<matrix_type> matrix( outer_length ) ;
  multi_vec_type x( Kokkos::ViewAllocateWithoutInitializing("x"),
                    inner_length, outer_length  ) ;
  multi_vec_type y( "y", inner_length, outer_length ) ;
  multi_vec_type tmp_x( "tmp_x", inner_length, outer_length ) ;
  multi_vec_type tmp_y( "tmp_y", inner_length, outer_length ) ;

  Kokkos::deep_copy( x , ScalarType(1.0) );

  for (size_t block=0; block<outer_length; ++block) {
    matrix_graph_type matrix_graph =
      Kokkos::create_staticcrsgraph<matrix_graph_type>(
        std::string("test crs graph") , fem_graph );

    matrix_values_type matrix_values = matrix_values_type(
      Kokkos::ViewAllocateWithoutInitializing("matrix"), graph_length );
    Kokkos::deep_copy(matrix_values , ScalarType(1.0) );
    matrix[block] = matrix_type("matrix", outer_length, matrix_values,
                                matrix_graph);
  }

  Device().fence();
  Kokkos::Timer clock ;
  int n_apply = 0;
  int n_add = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    typedef typename Cijk_type::k_iterator k_iterator;
    typedef typename Cijk_type::kj_iterator kj_iterator;
    typedef typename Cijk_type::kji_iterator kji_iterator;
    n_apply = 0;
    n_add = 0;
    k_iterator k_begin = Cijk->k_begin();
    k_iterator k_end = Cijk->k_end();
    for (k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      unsigned nj = Cijk->num_j(k_it);
      if (nj > 0) {
        int k = index(k_it);
        kj_iterator j_begin = Cijk->j_begin(k_it);
        kj_iterator j_end = Cijk->j_end(k_it);
        unsigned jdx = 0;
        for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          int j = index(j_it);
          vec_type xx = Kokkos::subview( x, Kokkos::ALL(), j );
          vec_type tt = Kokkos::subview( tmp_x, Kokkos::ALL(), jdx++ );
          Kokkos::deep_copy(tt, xx);
        }
        multi_vec_type tmp_x_view =
          Kokkos::subview( tmp_x, Kokkos::ALL(),
                                           std::make_pair(0u,nj));
        multi_vec_type tmp_y_view =
          Kokkos::subview( tmp_y, Kokkos::ALL(),
                                           std::make_pair(0u,nj));
        KokkosSparse::spmv( "N" , value_type(1.0) , matrix[k] , tmp_x_view , value_type(0.0) , tmp_y_view );
        n_apply += nj;
        jdx = 0;
        for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
          vec_type tmp_y_view =
            Kokkos::subview( tmp_y, Kokkos::ALL(), jdx++ );
          kji_iterator i_begin = Cijk->i_begin(j_it);
          kji_iterator i_end =  Cijk->i_end(j_it);
          for (kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
            int i = index(i_it);
            value_type c = value(i_it);
            vec_type y_view = Kokkos::subview( y, Kokkos::ALL(), i );
            //Stokhos::update( value_type(1.0) , y_view , c , tmp_y_view );
            KokkosBlas::update(value_type(1.0) , y_view, c, tmp_y_view, value_type(0.0), y_view);
            ++n_add;
          }
        }
      }
    }

  }
  Device().fence();

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
#endif

template< class Scalar, class Device >
void performance_test_driver_all( const int pdeg ,
                                  const int minvar ,
                                  const int maxvar ,
                                  const int nGrid ,
                                  const int nIter ,
                                  const bool test_block ,
                                  const bool symmetric )
{
  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_nonzeros =
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  std::cout.precision(8);

  //------------------------------

  std::cout << std::endl << "\"FEM NNZ = " << fem_nonzeros << "\"" << std::endl;

  std::cout << std::endl
            << "\"#nGrid\" , "
            << "\"#Variable\" , "
            << "\"PolyDegree\" , "
            << "\"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
            << "\"Original-Flat MXV-Time\" , "
            << "\"Original-Flat MXV-Speedup\" , "
            << "\"Original-Flat MXV-GFLOPS\" , "
            << "\"Commuted-Flat MXV-Speedup\" , "
            << "\"Commuted-Flat MXV-GFLOPS\" , "
            << "\"Block-Diagonal MXV-Speedup\" , "
            << "\"Block-Diagonal MXV-GFLOPS\" , "
            << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {

    std::vector<int> var_degree( nvar , pdeg );

    //------------------------------

    const std::vector<double> perf_flat_original =
      test_product_flat_original_matrix<Scalar,Device>(
        var_degree , nGrid , nIter , symmetric );

    const std::vector<double> perf_flat_commuted =
      test_product_flat_commuted_matrix<Scalar,Device>(
        var_degree , nGrid , nIter , symmetric );

    const std::vector<double> perf_matrix =
      test_product_tensor_diagonal_matrix<Scalar,Device>(
        var_degree , nGrid , nIter , symmetric );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(
        var_degree , nGrid , nIter , symmetric );

    if ( perf_flat_commuted[0] != perf_flat_original[0] ||
         perf_flat_commuted[3] != perf_flat_original[3] ) {
      std::cout << "ERROR: Original and commuted matrix sizes do not match"
                << std::endl
                << "  original size = " << perf_flat_original[0]
                << " , nonzero = " << perf_flat_original[3]
                << std::endl
                << "  commuted size = " << perf_flat_commuted[0]
                << " , nonzero = " << perf_flat_commuted[3]
                << std::endl ;
    }

    std::cout << nGrid << " , "
              << nvar << " , "
              << pdeg << " , "
              << perf_crs_tensor[4] << " , "
              << perf_crs_tensor[3] << " , "
              << perf_flat_original[0] << " , "
              << perf_flat_original[1] << " , "
              << perf_flat_original[1] / perf_flat_original[1] << " , "
              << perf_flat_original[2] << " , "
              << perf_flat_original[1] / perf_flat_commuted[1] << " , "
              << perf_flat_commuted[2] << " , "
              << perf_flat_original[1] / perf_matrix[1] << " , "
              << perf_matrix[2] << " , "
              << perf_flat_original[1] / perf_crs_tensor[1] << " , "
              << perf_crs_tensor[2] << " , "
              << std::endl ;
  }

  //------------------------------
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
            // << "\"Block-Coo-Tensor MXV-Speedup\" , "
            // << "\"Block-Coo-Tensor MXV-GFLOPS\" , "
            // << "\"Block-Tiled Crs-Tensor MXV-Speedup\" , "
            // << "\"Block-Tiled Crs-3-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Stokhos::CrsProductTensor<Scalar,Device>,Device>(
        var_degree , nGrid , nIter , symmetric );

    // const bool Pack = std::is_same<Device,Kokkos::Cuda>::value;
    // const std::vector<double> perf_coo_tensor =
    //   test_product_tensor_matrix<Scalar,Stokhos::CooProductTensor<Scalar,Device,Pack>,Device>(
    //     var_degree , nGrid , nIter , symmetric );

    // const std::vector<double> perf_tiled_crs_tensor =
    //   test_simple_tiled_product_tensor_matrix<Scalar,Device>(
    //     var_degree , nGrid , nIter , symmetric );

    std::vector<double> perf_original_mat_free_block;
#if defined(HAVE_STOKHOS_KOKKOSLINALG)
#if defined( KOKKOS_ENABLE_CUDA )
    enum { is_cuda = std::is_same<Device,Kokkos::Cuda>::value };
#else
    enum { is_cuda = false };
#endif
    if ( is_cuda )
      perf_original_mat_free_block =
        test_original_matrix_free_kokkos<Scalar,Device>(
          var_degree , nGrid , nIter , test_block , symmetric );
    else
      perf_original_mat_free_block =
        test_original_matrix_free_view<Scalar,Device,SparseMatOps>(
          var_degree , nGrid , nIter , test_block , symmetric );
#else
    perf_original_mat_free_block =
      test_original_matrix_free_view<Scalar,Device,SparseMatOps>(
        var_degree , nGrid , nIter , test_block , symmetric );
#endif

    std::cout << nGrid << " , "
              << nvar << " , "
              << pdeg << " , "
              << perf_crs_tensor[4] << " , "
              << perf_crs_tensor[3] << " , "
              << perf_original_mat_free_block[0] << " , "
              << perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[1] /
                 perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
              << perf_original_mat_free_block[1] / perf_crs_tensor[1] << " , "
              << perf_crs_tensor[2] << " , "
              // << perf_original_mat_free_block[1] / perf_coo_tensor[1] << " , "
              // << perf_coo_tensor[2] << " , "
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
    std::is_same<typename Device::memory_space,Kokkos::HostSpace>::value ;

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
struct performance_test_driver;

//----------------------------------------------------------------------------

}
