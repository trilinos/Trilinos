
#include <utility>
#include <cmath>
#include <iostream>

#include <KokkosArray_Host.hpp>
#include <KokkosArray_ProductTensor.hpp>
#include <KokkosArray_LegendrePolynomial.hpp>
#include <KokkosArray_SymmetricDiagonalSpec.hpp>
#include <KokkosArray_StochasticProductTensor.hpp>
#include <KokkosArray_CrsMatrix.hpp>
#include <KokkosArray_CrsProductTensorLegendre.hpp>
#include <KokkosArray_BlockCrsMatrix.hpp>

#include <impl/KokkosArray_Timer.hpp>

#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

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

template< typename ScalarType , class Device ,
          template< unsigned , typename , class > class TensorType >
std::vector<double>
test_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type** ,
                             KokkosArray::LayoutLeft ,
                             Device > block_vector_type ;

  typedef KokkosArray::NormalizedLegendrePolynomialBases<8> polynomial ;

  typedef KokkosArray::StochasticProductTensor< value_type , polynomial , Device , TensorType > tensor_type ;

  typedef KokkosArray::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  matrix.block = KokkosArray::create_product_tensor< tensor_type >( var_degree );
  matrix.graph = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();

  matrix.values = block_vector_type( "matrix" , inner_length , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM = 
    KokkosArray::create_mirror( matrix.values );

  for ( size_t i=0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1.0;
    }
  }
  
  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:
  
  typename block_vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1.0 ;
    }
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = 5.0 * matrix.block.tensor().num_non_zeros() +
    inner_length;
  const double flops = 1.0e-9*graph_length*flops_per_block / seconds_per_iter;

  std::vector<double> perf(5) ;

  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops;
  perf[3] = matrix.block.tensor().entry_count();
  perf[4] = inner_length ;

  return perf ;
}

//----------------------------------------------------------------------------

template< typename VectorScalar , typename MatrixScalar , typename TensorScalar , class Device >
std::vector<double>
test_product_tensor_legendre(
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount )
{
  typedef KokkosArray::View< VectorScalar** ,
                             KokkosArray::LayoutLeft ,
                             Device > vector_type ;

  typedef KokkosArray::CrsProductTensorLegendre< TensorScalar , Device >  tensor_type ;

  typedef KokkosArray::BlockCrsMatrix< tensor_type , MatrixScalar , Device > matrix_type ;

  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;
  const size_t fem_graph_length = unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  const std::vector<unsigned> var_degree( arg_var_degree.begin() , arg_var_degree.end() );

  const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
    tensor( var_degree );

  const size_t stoch_length = tensor.bases_count();

  std::vector< std::vector< size_t > > stoch_graph( stoch_length );

  for ( size_t i = 0 ; i < stoch_length ; ++i ) {
    for ( size_t j = 0 ; j < stoch_length ; ++j ) {
      if ( KokkosArray::matrix_nonzero(tensor,i,j) ) {
        stoch_graph[i].push_back(j);
      }
    }
  }

  //------------------------------
  // Generate input multivector:
  
  vector_type x = vector_type( "x" , stoch_length , fem_length );
  vector_type y = vector_type( "y" , stoch_length , fem_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  
  for ( size_t iColFEM = 0 ;   iColFEM < fem_length ;   ++iColFEM ) {
  for ( size_t iColStoch = 0 ; iColStoch < stoch_length ; ++iColStoch ) {
    hx(iColStoch,iColFEM) = 1.0;
  }}

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.block = tensor_type( var_degree );

  matrix.graph = KokkosArray::create_crsarray<graph_type>( std::string("test crs graph") , fem_graph );

  if ( stoch_length != matrix.block.dimension() ) {
    throw std::runtime_error("test_product_tensor_legendre matrix sizing error");
  }

  matrix.values = vector_type( "matrix" , stoch_length , fem_graph_length );

  typename vector_type::HostMirror hM = KokkosArray::create_mirror( matrix.values );

  for ( size_t iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {
    for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
      
      for ( size_t k = 0 ; k < stoch_length ; ++k ) {
        hM(k,iEntryFEM) = 1.0;
      }
    }
  }

  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------

  const KokkosArray::Impl::Multiply< matrix_type , vector_type , vector_type > op( matrix , x , y );

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    op.run();
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops_per_block = matrix.block.multiply_add_flops();
  const double flops = 1.0e-9*fem_graph_length*flops_per_block / seconds_per_iter;

  std::vector<double> perf(3) ;

  perf[0] = fem_length * stoch_length ;
  perf[1] = seconds_per_iter ;
  perf[2] = flops ;

  return perf ;
}

//----------------------------------------------------------------------------

template< typename ScalarType , class Device >
std::vector<double>
test_product_tensor_diagonal_matrix(
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type**,
                             KokkosArray::LayoutLeft ,
                             Device > block_vector_type ;

  //------------------------------

  typedef KokkosArray::BlockCrsMatrix< KokkosArray::SymmetricDiagonalSpec< Device > ,
                                  value_type , Device > matrix_type ;

  typedef typename matrix_type::graph_type  graph_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;
  const size_t fem_graph_length = unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const std::vector<unsigned> var_degree( arg_var_degree.begin() , arg_var_degree.end() );

  const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
    tensor( var_degree );

  const size_t stoch_length = tensor.bases_count();

  //------------------------------

  block_vector_type x = block_vector_type( "x" , stoch_length , fem_length );
  block_vector_type y = block_vector_type( "y" , stoch_length , fem_length );

  typename block_vector_type::HostMirror hx        = KokkosArray::create_mirror( x );

  for ( size_t iColFEM = 0 ;   iColFEM < fem_length ;   ++iColFEM ) {
  for ( size_t iColStoch = 0 ; iColStoch < stoch_length ; ++iColStoch ) {
    hx(iColStoch,iColFEM) = 1.0;
  }}

  KokkosArray::deep_copy( x , hx );

  //------------------------------
  // Generate CRS matrix of blocks with symmetric diagonal storage

  matrix_type matrix ;

  matrix.block  = KokkosArray::SymmetricDiagonalSpec< Device >( stoch_length );
  matrix.graph  = KokkosArray::create_crsarray<graph_type>( std::string("test product tensor graph") , fem_graph );
  matrix.values = block_vector_type( "matrix" , matrix.block.matrix_size() , fem_graph_length );

  {
    typename block_vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRowStoch = 0 ; iRowStoch < stoch_length ; ++iRowStoch ) {
      for ( size_t iRowFEM = 0 , iEntryFEM = 0 ; iRowFEM < fem_length ; ++iRowFEM ) {

        for ( size_t iRowEntryFEM = 0 ; iRowEntryFEM < fem_graph[iRowFEM].size() ; ++iRowEntryFEM , ++iEntryFEM ) {
          
          for ( size_t iColStoch = 0 ; iColStoch < stoch_length ; ++iColStoch ) {

            const size_t offset = matrix.block.matrix_offset( iRowStoch , iColStoch );

            hM( offset , iEntryFEM ) = 1.0 ;
          }
        }

      }
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*fem_graph_length*stoch_length*stoch_length / seconds_per_iter;

  std::vector<double> perf(3);
  perf[0] = fem_length * stoch_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
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
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type[] , Device > vector_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;

  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Tensor evaluation:

  const std::vector<unsigned> var_degree( arg_var_degree.begin() , arg_var_degree.end() );

  const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
    tensor( var_degree );

  const size_t stoch_length = tensor.bases_count();

  std::vector< std::vector< size_t > > stoch_graph( stoch_length );

  for ( size_t i = 0 ; i < stoch_length ; ++i ) {
    for ( size_t j = 0 ; j < stoch_length ; ++j ) {
      if ( KokkosArray::matrix_nonzero(tensor,i,j) ) {
        stoch_graph[i].push_back(j);
      }
    }
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

  vector_type x = vector_type( "x" , flat_length );
  vector_type y = vector_type( "y" , flat_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  
  for ( size_t iCol = 0 ; iCol < flat_length ; ++iCol ) {
    hx(iCol) = 1.0;
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values = vector_type( "matrix" , flat_graph_length );
  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRow = 0 , iEntry = 0 ; iRow < flat_length ; ++iRow ) {
      
      for ( size_t iRowEntry = 0 ; iRowEntry < flat_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {

        hM( iEntry ) = 1.0 ;
      }

    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //KokkosArray::write_matrix_market(matrix, "flat_commuted.mm");

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

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
  const std::vector<int> & arg_var_degree ,
  const int nGrid ,
  const int iterCount )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type[] , Device > vector_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t fem_length = nGrid * nGrid * nGrid ;

  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Tensor evaluation:

  const std::vector<unsigned> var_degree( arg_var_degree.begin() , arg_var_degree.end() );

  const KokkosArray::TripleProductTensorLegendreCombinatorialEvaluation
    tensor( var_degree );

  const size_t stoch_length = tensor.bases_count();

  std::vector< std::vector< size_t > > stoch_graph( stoch_length );

  for ( size_t i = 0 ; i < stoch_length ; ++i ) {
    for ( size_t j = 0 ; j < stoch_length ; ++j ) {
      if ( KokkosArray::matrix_nonzero(tensor,i,j) ) {
        stoch_graph[i].push_back(j);
      }
    }
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

  vector_type x = vector_type( "x" , flat_length );
  vector_type y = vector_type( "y" , flat_length );

  typename vector_type::HostMirror hx        = KokkosArray::create_mirror( x );
  
  for ( size_t iCol = 0 ; iCol < flat_length ; ++iCol ) {
    hx(iCol) = 1.0;
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values = vector_type( "matrix" , flat_graph_length );
  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t iRow = 0 , iEntry = 0 ; iRow < flat_length ; ++iRow ) {

      for ( size_t iRowEntry = 0 ; iRowEntry < flat_graph[ iRow ].size() ; ++iRowEntry , ++iEntry ) {

        hM( iEntry ) = 1.0 ;

      }

    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //KokkosArray::write_matrix_market(matrix, "flat_original.mm");

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

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
test_original_matrix_free_vec(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool test_block )
{
  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,value_type> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++)
    bases[i] = Teuchos::rcp(new basis_type(var_degree[i],true));
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(
	  bases, ScalarTolerances<value_type>::sparse_cijk_tol()));
  const size_t outer_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor(outer_length);

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = 
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  
  typedef KokkosArray::View<value_type[],Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  std::vector<vec_type> x( outer_length ) ;
  std::vector<vec_type> y( outer_length ) ;
  std::vector<vec_type> tmp( outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

    matrix[block].values = vec_type( "matrix" , graph_length );

    x[block]   = vec_type( "x" , inner_length );
    y[block]   = vec_type( "y" , inner_length );
    tmp[block] = vec_type( "tmp" , inner_length );

    typename vec_type::HostMirror hM =
      KokkosArray::create_mirror( matrix[block].values );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1.0 ;
    }

    KokkosArray::deep_copy( matrix[block].values , hM );

    typename vec_type::HostMirror hx =
      KokkosArray::create_mirror( x[block] );

    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      hx(i) = 1.0 ;
    }

    KokkosArray::deep_copy( x[block] , hx );
  }
  

  KokkosArray::Impl::Timer clock ;
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
        KokkosArray::multiply( matrix[k] , xx , yy, test_block );
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
	    KokkosArray::update( value_type(1.0) , y[i] , c , yy[jdx] );
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

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length;
  perf[1] = seconds_per_iter ;
  perf[2] = flops/seconds_per_iter;
  perf[3] = flops;

  return perf;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class Scalar, class Device >
void performance_test_driver_all( const int pdeg ,
				  const int minvar ,
				  const int maxvar ,
				  const int nGrid ,
				  const int nIter ,
				  const bool test_block )
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
	    << "\"Block-Legendre-Tensor MXV-Speedup\" , "
            << "\"Block-Legendre-Tensor MXV-GFLOPS\" , "
	    << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {

    std::vector<int> var_degree( nvar , pdeg );

    //------------------------------

    const std::vector<double> perf_flat_original =
      test_product_flat_original_matrix<Scalar,Device>( 
	var_degree , nGrid , nIter );

    const std::vector<double> perf_flat_commuted =
      test_product_flat_commuted_matrix<Scalar,Device>( 
	var_degree , nGrid , nIter );

    const std::vector<double> perf_matrix =
      test_product_tensor_diagonal_matrix<Scalar,Device>( 
	var_degree , nGrid , nIter );

    const std::vector<double> perf_legendre_tensor =
      test_product_tensor_legendre<Scalar,Scalar,Scalar,Device>( 
	var_degree , nGrid , nIter );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Device,KokkosArray::CrsProductTensor>( 
	var_degree , nGrid , nIter );

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
	      << perf_flat_original[1] / perf_legendre_tensor[1] << " , "
	      << perf_legendre_tensor[2] << " , "
	      << perf_flat_original[1] / perf_crs_tensor[1] << " , "
	      << perf_crs_tensor[2] << " , "
	      << std::endl ;
  }

  //------------------------------
}

template< class Scalar, class Device >
void performance_test_driver_poly( const int pdeg ,
				   const int minvar ,
				   const int maxvar ,
				   const int nGrid ,
				   const int nIter ,
				   const bool test_block )
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
	    << "\"Block-Legendre-Tensor MXV-Speedup\" , "
            << "\"Block-Legendre-Tensor MXV-GFLOPS\" , "
	    << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const std::vector<double> perf_crs_tensor =
      test_product_tensor_matrix<Scalar,Device,KokkosArray::CrsProductTensor>( 
	var_degree , nGrid , nIter );

    const std::vector<double> perf_legendre_tensor =
      test_product_tensor_legendre<Scalar,Scalar,Scalar,Device>(
	var_degree , nGrid , nIter );

    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<Scalar,Device>( 
	var_degree , nGrid , nIter , test_block );

    std::cout << nGrid << " , "
	      << nvar << " , " 
	      << pdeg << " , "
	      << perf_crs_tensor[4] << " , "
	      << perf_crs_tensor[3] << " , "
	      << perf_original_mat_free_block[0] << " , "
	      << perf_original_mat_free_block[1] << " , "
	      << perf_original_mat_free_block[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
	      << perf_original_mat_free_block[1] / perf_legendre_tensor[1] << " , "
              << perf_legendre_tensor[2] << " , "
	      << perf_original_mat_free_block[1] / perf_crs_tensor[1] << " , "
              << perf_crs_tensor[2] << " , "

	      << std::endl ;
  }

  //------------------------------
}

template< class Scalar, class Device >
struct performance_test_driver {
  static void run(bool test_flat, bool test_orig, bool test_block) {}
};

//----------------------------------------------------------------------------

}


