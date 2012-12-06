
#include <utility>
#include <cmath>
#include <iostream>

#include <impl/KokkosArray_Timer.hpp>

#ifdef HAVE_KOKKOSARRAY_STOKHOS
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#endif

namespace unit_test {

template< unsigned P >
void test_integration()
{
  KokkosArray::GaussLegendre<P> rule ;
  double result_1 = 0 ;
  double result_x = 0 ;
  double result_x2 = 0 ;
  for ( unsigned i = 0 ; i < rule.N ; ++i ) {
    result_1 += rule.weights[i];
    result_x += rule.points[i] * rule.weights[i];
    result_x2 += rule.points[i] * rule.points[i] * rule.weights[i];
  }
  std::cout << "IntegrateP" << P << "(1) = " << result_1 << std::endl ;
  std::cout << "IntegrateP" << P << "(x) = " << result_x << std::endl ;
  std::cout << "IntegrateP" << P << "(x^2) = " << result_x2 << std::endl ;
}

//----------------------------------------------------------------------------

template< unsigned P , class Device >
void test_inner_product_legengre_polynomial()
{
  const double tolerance = 1e-14 ;

  KokkosArray::GaussLegendre<P*2> rule ;
  KokkosArray::NormalizedLegendrePolynomialBases<P> poly ;

  double values[P+1];
  double result[P+1][P+1];

  for ( unsigned k = 0 ; k <= P ; ++k ) {
    for ( unsigned j = 0 ; j <= P ; ++j ) {
      result[k][j] = 0 ;
    }
  }

  for ( unsigned i = 0 ; i < rule.N ; ++i ) {
    poly.evaluate( P , rule.points[i] , values );

    for ( unsigned k = 0 ; k <= P ; ++k ) {
      for ( unsigned j = 0 ; j <= P ; ++j ) {
        result[k][j] += rule.weights[i] * values[k] * values[j] ;
      }
    }
  }

  for ( unsigned k = 0 ; k <= P ; ++k ) {
    if ( tolerance < std::fabs( result[k][k] ) ) {
      std::cout << "<P" << k << ",P" << k << "> = " ;
      std::cout.precision(16);
      std::cout << result[k][k] << std::endl ;
    }
  }

  for ( unsigned k = 0 ; k <= P ; ++k ) {
    for ( unsigned j = k + 1 ; j <= P ; ++j ) {
      if ( tolerance < std::fabs( result[k][j] ) ) {
        std::cout << "<P" << k << ",P" << j << "> = " ;
        std::cout.precision(16);
        std::cout << result[k][j] << std::endl ;
      }
    }
  }
}

//----------------------------------------------------------------------------

template< unsigned P , class Device >
void test_triple_product_legendre_polynomial()
{
  const double tolerance = 1e-14 ;

  KokkosArray::GaussLegendre<P*3+1> rule ;
  KokkosArray::NormalizedLegendrePolynomialBases<P> poly ;

  double values[P+1];
  double result[P+1][P+1][P+1];

  for ( unsigned k = 0 ; k <= P ; ++k ) {
  for ( unsigned j = 0 ; j <= P ; ++j ) {
  for ( unsigned i = 0 ; i <= P ; ++i ) {
      result[k][j][i] = 0 ;
  } } }

  for ( unsigned n = 0 ; n < rule.N ; ++n ) {
    poly.evaluate( P , rule.points[n] , values );

    for ( unsigned k = 0 ; k <= P ; ++k ) {
    for ( unsigned j = 0 ; j <= P ; ++j ) {
    for ( unsigned i = 0 ; i <= P ; ++i ) {
      result[k][j][i] += rule.weights[n] * values[k] * values[j] * values[i] ;
    } } }
  }

  for ( unsigned k = 0 ; k <= P ; ++k ) {
    if ( tolerance < std::fabs( result[k][k][k] ) ) {
      std::cout << "<P" << k << ",P" << k << ",P" << k << "> = " ;
      std::cout.precision(16);
      std::cout << result[k][k][k] << std::endl ;
    }
  }

  for ( unsigned k = 0 ; k <= P ; ++k ) {
  for ( unsigned j = k + 1 ; j <= P ; ++j ) {
  for ( unsigned i = j ; i <= P ; ++i ) {
    if ( tolerance < std::fabs( result[k][j][i] ) ) {
      std::cout << "<P" << k << ",P" << j << ",P" << i << "> = " ;
      std::cout.precision(16);
      std::cout << result[k][j][i] << std::endl ;
    }
  } } }
}

//----------------------------------------------------------------------------

template< class Device , template< unsigned , typename , class > class TensorType >
void test_product_tensor( const std::vector<int> & var_degree )
{
  typedef KokkosArray::NormalizedLegendrePolynomialBases<4> polynomial ;
  typedef KokkosArray::StochasticProductTensor< double , polynomial , Device , TensorType > tensor_type ;

  tensor_type tensor = KokkosArray::create_product_tensor< tensor_type >( var_degree );

  // Verification?
}

//----------------------------------------------------------------------------

template< typename ScalarType , class Device ,
          template< unsigned , typename , class > class TensorType >
std::pair<size_t,double>
test_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type** ,
                             KokkosArray::LayoutLeft ,
                             Device > block_vector_type ;

  typedef KokkosArray::NormalizedLegendrePolynomialBases<4> polynomial ;

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
  const size_t inner_matrix_size = matrix.block.dimension();

  matrix.values = block_vector_type( "matrix" , inner_matrix_size , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM = KokkosArray::create_mirror( matrix.values );
  
  for ( size_t i = 0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1 + i ;
    }
  }
  
  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:
  
  typename block_vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
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

  //------------------------------

  if ( print_flag ) {
    typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );

    KokkosArray::deep_copy( hy , y );

    std::cout << std::endl << "test_product_tensor_matrix" << std::endl ;
    for ( size_t i = 0 ; i < outer_length ; ++i ) {
      std::cout << "hy(:," << i << ") =" ;
      for ( size_t j = 0 ; j < inner_length ; ++j ) {
        std::cout << " " << hy(j,i);
      }
      std::cout << std::endl ;
    }
  }

  return std::pair<size_t,double>( outer_length * inner_length , seconds_per_iter );
}

//----------------------------------------------------------------------------

template< typename ScalarType , class Device >
std::vector<double>
test_product_tensor_diagonal_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type**,
                             KokkosArray::LayoutLeft ,
                             Device > block_vector_type ;

  typedef KokkosArray::NormalizedLegendrePolynomialBases<4> polynomial ;
  typedef KokkosArray::StochasticProductTensor< value_type , polynomial , KokkosArray::Host , KokkosArray::SparseProductTensor > tensor_type ;

  //------------------------------

  typedef KokkosArray::Impl::Multiply<
            typename tensor_type::tensor_type ,
            KokkosArray::SymmetricDiagonalSpec< KokkosArray::Host > ,
            void > multiply_type ;

  typedef KokkosArray::BlockCrsMatrix< KokkosArray::SymmetricDiagonalSpec< Device > ,
                                  value_type , Device > matrix_type ;

  typedef typename matrix_type::graph_type  graph_type ;
  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const tensor_type tensor =
    KokkosArray::create_product_tensor< tensor_type >( var_degree );

  const size_t inner_length = tensor.dimension();

  //------------------------------
  // Generate CRS matrix of blocks with symmetric diagonal storage

  matrix_type matrix ;

  matrix.block  = KokkosArray::SymmetricDiagonalSpec< Device >( inner_length );
  matrix.graph  = KokkosArray::create_crsarray<graph_type>( std::string("test product tensor graph") , graph );
  matrix.values = block_vector_type( "matrix" , matrix.block.matrix_size() , graph_length );

  block_vector_type x = block_vector_type( "x" , inner_length , outer_length );
  block_vector_type y = block_vector_type( "y" , inner_length , outer_length );

  typename block_vector_type::HostMirror hM =
    KokkosArray::create_mirror( matrix.values );

  {
    std::vector< value_type > a( inner_length );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      for ( size_t j = 0 ; j < inner_length ; ++j ) {
        a[j] = 1 + j + 10 * i ;
      }
      // Tensor expansion:
      multiply_type::apply( tensor.tensor() , & a[0] , & hM(0,i) );
    }
  }

  KokkosArray::deep_copy( matrix.values , hM );

  //------------------------------

  typename block_vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
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
  const double flops = 2.0*1e-9*graph_length*inner_length*inner_length / seconds_per_iter;

  //------------------------------

  if ( print_flag ) {
    typename block_vector_type::HostMirror hy = KokkosArray::create_mirror( y );

    KokkosArray::deep_copy( hy , y );

    std::cout << std::endl << "test_product_tensor_diagonal_matrix"
              << std::endl ;
    for ( size_t i = 0 ; i < outer_length ; ++i ) {
      std::cout << "hy(:," << i << ") =" ;
      for ( size_t j = 0 ; j < inner_length ; ++j ) {
        std::cout << " " << hy(j,i);
      }
      std::cout << std::endl ;
    }
  }

  std::vector<double> perf(3);
  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  return perf;
}

//----------------------------------------------------------------------------
// Flatten to a plain CRS matrix

template< typename ScalarType , class Device >
std::vector<double>
test_product_flat_commuted_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type[] , Device > vector_type ;

  typedef KokkosArray::NormalizedLegendrePolynomialBases<4> polynomial ;

  typedef KokkosArray::StochasticProductTensor<
     value_type , polynomial ,
     KokkosArray::Host ,
     KokkosArray::CrsProductTensor > tensor_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const tensor_type tensor =
    KokkosArray::create_product_tensor< tensor_type >( var_degree );

  const size_t inner_length = tensor.dimension();

  std::vector< std::vector<size_t> > tensor_graph( inner_length );

  for ( size_t iInnerRow = 0 ; iInnerRow < inner_length ; ++iInnerRow ) {
  for ( size_t iInnerCol = 0 ; iInnerCol < inner_length ; ++iInnerCol ) {

    for ( KokkosArray::Host::size_type
            n = tensor.tensor().entry_begin( iInnerRow ) ;
            n < tensor.tensor().entry_end(   iInnerRow ) ; ++n ) {

      if ( iInnerCol == tensor.tensor().coord( n , 0 ) ||
           iInnerCol == tensor.tensor().coord( n , 1 ) ) {
        tensor_graph[iInnerRow].push_back( iInnerCol );
        break ;
      }
    }
  }
  }

  //------------------------------
  // Generate flattened graph:
  //
  // dof(i,j) -> dof(i+j*inner_length)

  const size_t flat_length  = inner_length * outer_length ;

  std::vector< std::vector<size_t> > flat_graph( flat_length );

  for ( size_t iOuterRow = 0 ; iOuterRow < outer_length ; ++iOuterRow ) {
  for ( size_t iInnerRow = 0 ; iInnerRow < inner_length ; ++iInnerRow ) {

    const size_t iFlatRow = iInnerRow + iOuterRow * inner_length ;

    const size_t iOuterNZ = fem_graph[iOuterRow].size();
    const size_t iInnerNZ = tensor_graph[iInnerRow].size();
    const size_t iFlatNZ  = iOuterNZ * iInnerNZ ;

    flat_graph[iFlatRow].resize( iFlatNZ );

    for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterNZ ; ++iOuterEntry ) {
    for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerNZ ; ++iInnerEntry ) {

      const size_t iOuterCol = fem_graph[   iOuterRow][iOuterEntry];
      const size_t iInnerCol = tensor_graph[iInnerRow][iInnerEntry];

      const size_t iFlatColumn = iInnerCol +   iOuterCol   * inner_length ;
      const size_t iFlatEntry  = iInnerEntry + iOuterEntry * iInnerNZ ;

      flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;
    }
    }
  }
  }

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values = vector_type( "matrix" , flat_graph_length );
  vector_type x = vector_type( "x" , flat_length );
  vector_type y = vector_type( "y" , flat_length );

  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t i = 0 ; i < flat_graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //------------------------------

  typename vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < flat_length ; ++i ) {
    hx(i) = 1 + i ;
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*flat_graph_length / seconds_per_iter;

  //------------------------------

  if ( print_flag ) {
    typename vector_type::HostMirror hy = KokkosArray::create_mirror( y );

    KokkosArray::deep_copy( hy , y );

    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
    for ( size_t i = 0 ; i < flat_length ; ++i ) {
      std::cout << "hy(" << i << ") = " << hy(i) << std::endl ;
    }
  }

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  perf[3] = flat_graph_length;
  return perf;
}

//----------------------------------------------------------------------------
// Flatten to a plain CRS matrix

template< typename ScalarType , class Device >
std::vector<double>
test_product_flat_original_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type[] , Device > vector_type ;

  typedef KokkosArray::NormalizedLegendrePolynomialBases<4> polynomial ;

  typedef KokkosArray::StochasticProductTensor<
    value_type , polynomial ,
    KokkosArray::Host ,
    KokkosArray::CrsProductTensor > tensor_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const tensor_type tensor =
    KokkosArray::create_product_tensor< tensor_type >( var_degree );

  const size_t inner_length = tensor.dimension();

  std::vector< std::vector<size_t> > tensor_graph( inner_length );

  for ( size_t iInnerRow = 0 ; iInnerRow < inner_length ; ++iInnerRow ) {
  for ( size_t iInnerCol = 0 ; iInnerCol < inner_length ; ++iInnerCol ) {

    for ( KokkosArray::Host::size_type
            n = tensor.tensor().entry_begin( iInnerRow ) ;
            n < tensor.tensor().entry_end(   iInnerRow ) ; ++n ) {

      if ( iInnerCol == tensor.tensor().coord( n , 0 ) ||
           iInnerCol == tensor.tensor().coord( n , 1 ) ) {
        tensor_graph[iInnerRow].push_back( iInnerCol );
        break ;
      }
    }
  }
  }

  //------------------------------
  // Generate flattened graph:
  //
  // dof(i,j) -> dof(i+j*inner_length)

  const size_t flat_length  = inner_length * outer_length ;

  std::vector< std::vector<size_t> > flat_graph( flat_length );

  for ( size_t iInnerRow = 0 ; iInnerRow < inner_length ; ++iInnerRow ) {
  for ( size_t iOuterRow = 0 ; iOuterRow < outer_length ; ++iOuterRow ) {

    const size_t iFlatRow = iOuterRow + iInnerRow * outer_length ;

    const size_t iOuterNZ = fem_graph[iOuterRow].size();
    const size_t iInnerNZ = tensor_graph[iInnerRow].size();
    const size_t iFlatNZ  = iOuterNZ * iInnerNZ ;

    flat_graph[iFlatRow].resize( iFlatNZ );

    for ( size_t iInnerEntry = 0 ; iInnerEntry < iInnerNZ ; ++iInnerEntry ) {
    for ( size_t iOuterEntry = 0 ; iOuterEntry < iOuterNZ ; ++iOuterEntry ) {

      const size_t iOuterCol = fem_graph[   iOuterRow][iOuterEntry];
      const size_t iInnerCol = tensor_graph[iInnerRow][iInnerEntry];

      const size_t iFlatColumn = iOuterCol   + iInnerCol   * outer_length ;
      const size_t iFlatEntry  = iOuterEntry + iInnerEntry * iOuterNZ ;

      flat_graph[iFlatRow][iFlatEntry] = iFlatColumn ;
    }
    }
  }
  }

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values = vector_type( "matrix" , flat_graph_length );

  vector_type x = vector_type( "x" , flat_length );
  vector_type y = vector_type( "y" , flat_length );

  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t i = 0 ; i < flat_graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //KokkosArray::write_matrix_market(matrix, "flat_original.mm");

  //------------------------------

  typename vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < flat_length ; ++i ) {
    hx(i) = 1 + i ;
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1e-9*flat_graph_length / seconds_per_iter;

  //------------------------------

  if ( print_flag ) {
    typename vector_type::HostMirror hy = KokkosArray::create_mirror( y );

    KokkosArray::deep_copy( hy , y );

    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
    for ( size_t i = 0 ; i < flat_length ; ++i ) {
      std::cout << "hy(" << i << ") = " << hy(i) << std::endl ;
    }
  }

  std::vector<double> perf(3);
  perf[0] = outer_length * inner_length ;
  perf[1] = seconds_per_iter;
  perf[2] = flops;
  return perf;
}

//----------------------------------------------------------------------------
// Outer original matrix-free block algorithm
#ifdef HAVE_KOKKOSARRAY_STOKHOS
template< typename ScalarType , class Device >
std::vector<double>
test_original_matrix_free_block(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false ,
  const bool test_block = false )
{
  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++)
    bases[i] = Teuchos::rcp(new basis_type(var_degree[i],true));
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(bases, 1e-12));
  const size_t outer_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor(outer_length);

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  
  typedef KokkosArray::View<value_type**, KokkosArray::LayoutLeft, Device> multi_vec_type ;
  typedef KokkosArray::View<value_type[], Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  multi_vec_type x( "x" , inner_length , outer_length ) ;
  multi_vec_type y( "y" , inner_length , outer_length ) ;
  multi_vec_type tmp( "tmp" , inner_length ,  outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );
    const size_t graph_length = matrix[block].graph.entries.dimension(0);
    matrix[block].values = vec_type( "matrix" , graph_length );

    typename vec_type::HostMirror hM =
      KokkosArray::create_mirror( matrix[block].values );
    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }
    KokkosArray::deep_copy( matrix[block].values , hM );
  }
  
  typename multi_vec_type::HostMirror hx =
    KokkosArray::create_mirror( x );
  for (size_t block=0; block<outer_length; ++block) {
    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      hx( i , block) = 1 + i ;
    }
  }
  KokkosArray::deep_copy( x , hx );

  KokkosArray::Impl::Timer clock ;
  int n_apply = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    n_apply = 0;
    Cijk_type::k_iterator k_begin = Cijk->k_begin();
    Cijk_type::k_iterator k_end = Cijk->k_end();
    for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int nj = Cijk->num_j(k_it);
      if (nj > 0) {
	int k = index(k_it);
	Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
	Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
	std::vector<int> j_indices(nj);
	int jdx = 0;
	for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	  j_indices[jdx++] = index(j_it);
	}
        KokkosArray::multiply( matrix[k] , x , tmp, j_indices , test_block );
        n_apply += nj;
	jdx = 0;
	for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	  const vec_type tmp_view( tmp , jdx++ );
	  Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
	  Cijk_type::kji_iterator i_end =  Cijk->i_end(j_it);
	  for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
	    int i = index(i_it);
	    value_type c = value(i_it);
	    const vec_type y_view( y , i );
	    KokkosArray::update( 1.0 , y_view , c , tmp_view );
	  }
	}
      }
    }

  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1.0e-9*n_apply*matrix[0].graph.entries.dimension(0);

  //------------------------------

  if ( print_flag ) {
    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
    typename multi_vec_type::HostMirror hy = KokkosArray::create_mirror( y );
    KokkosArray::deep_copy( hy , y );
    std::cout << "hy = " << std::endl ;
    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      for ( size_t j = 0 ; j < outer_length ; ++j ) {
	std::cout << " " << hy(i,j);
      }
      std::cout << std::endl ;
    }
  }

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length;
  perf[1] = seconds_per_iter ;
  perf[2] = flops/seconds_per_iter;
  perf[3] = flops;

  return perf;
}

template< typename ScalarType , class Device >
std::vector<double>
test_original_matrix_free_vec(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false ,
  const bool test_block = false )
{
  typedef ScalarType value_type ;
  typedef Stokhos::OneDOrthogPolyBasis<int,value_type> abstract_basis_type;
  typedef Stokhos::LegendreBasis<int,value_type> basis_type;
  typedef Stokhos::CompletePolynomialBasis<int,value_type> product_basis_type;
  typedef Stokhos::Sparse3Tensor<int,double> Cijk_type;

  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Create Stochastic Galerkin basis and expansion
  const size_t num_KL = var_degree.size();
  Array< RCP<const abstract_basis_type> > bases(num_KL); 
  for (size_t i=0; i<num_KL; i++)
    bases[i] = Teuchos::rcp(new basis_type(var_degree[i],true));
  RCP<const product_basis_type> basis = 
    rcp(new product_basis_type(bases, 1e-12));
  const size_t outer_length = basis->size();
  RCP<Cijk_type> Cijk = basis->computeTripleProductTensor(outer_length);

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t inner_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  
  typedef KokkosArray::View<value_type[],Device> vec_type ;

  std::vector<matrix_type> matrix( outer_length ) ;
  std::vector<vec_type> x( outer_length ) ;
  std::vector<vec_type> y( outer_length ) ;
  std::vector<vec_type> tmp( outer_length ) ;

  for (size_t block=0; block<outer_length; ++block) {
    matrix[block].graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

    const size_t graph_length = matrix[block].graph.entries.dimension(0);

    matrix[block].values = vec_type( "matrix" , graph_length );

    x[block]   = vec_type( "x" , inner_length );
    y[block]   = vec_type( "y" , inner_length );
    tmp[block] = vec_type( "tmp" , inner_length );

    typename vec_type::HostMirror hM =
      KokkosArray::create_mirror( matrix[block].values );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    KokkosArray::deep_copy( matrix[block].values , hM );

    typename vec_type::HostMirror hx =
      KokkosArray::create_mirror( x[block] );

    for ( size_t i = 0 ; i < inner_length ; ++i ) {
      hx(i) = 1 + i ;
    }

    KokkosArray::deep_copy( x[block] , hx );
  }
  

  KokkosArray::Impl::Timer clock ;
  int n_apply = 0;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {

    // Original matrix-free multiply algorithm using a block apply
    n_apply = 0;
    Cijk_type::k_iterator k_begin = Cijk->k_begin();
    Cijk_type::k_iterator k_end = Cijk->k_end();
    for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      int nj = Cijk->num_j(k_it);
      if (nj > 0) {
	int k = index(k_it);
	Cijk_type::kj_iterator j_begin = Cijk->j_begin(k_it);
	Cijk_type::kj_iterator j_end = Cijk->j_end(k_it);
	std::vector<vec_type> xx(nj), yy(nj);
	int jdx = 0;
	for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	  int j = index(j_it);
	  xx[jdx] = x[j];
	  yy[jdx] = tmp[j];
	  jdx++;
	}
        KokkosArray::multiply( matrix[k] , xx , yy, test_block );
        n_apply += nj;
	jdx = 0;
	for (Cijk_type::kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
	  Cijk_type::kji_iterator i_begin = Cijk->i_begin(j_it);
	  Cijk_type::kji_iterator i_end =  Cijk->i_end(j_it);
	  for (Cijk_type::kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
	    int i = index(i_it);
	    value_type c = value(i_it);
	    KokkosArray::update( 1.0 , y[i] , c , yy[jdx] );
	  }
	  jdx++;
	}
      }
    }

  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );
  const double flops = 2.0*1.0e-9*n_apply*matrix[0].graph.entries.dimension(0);

  //------------------------------

  if ( print_flag ) {
    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
    for ( size_t i = 0 ; i < outer_length ; ++i ) {
      typename vec_type::HostMirror hy = KokkosArray::create_mirror( y[i] );

      KokkosArray::deep_copy( hy , y[i] );

      std::cout << "hy(:," << i << ") =" ;
      for ( size_t j = 0 ; j < inner_length ; ++j ) {
        std::cout << " " << hy(j);
      }
      std::cout << std::endl ;
    }
  }

  std::vector<double> perf(4);
  perf[0] = outer_length * inner_length;
  perf[1] = seconds_per_iter ;
  perf[2] = flops/seconds_per_iter;
  perf[3] = flops;

  return perf;
}
#endif

//----------------------------------------------------------------------------
// A plain CRS matrix

template< typename ScalarType , class Device >
std::pair<size_t,double>
test_flat_matrix(
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef KokkosArray::View< value_type[] , Device > vector_type ;

  //------------------------------

  typedef KokkosArray::CrsMatrix<value_type,Device> matrix_type ;
  typedef KokkosArray::CrsArray<int,Device,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t length = nGrid * nGrid * nGrid ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = KokkosArray::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

  matrix.values = vector_type( "matrix" , graph_length );
  vector_type x = vector_type( "x" , length );
  vector_type y = vector_type( "y" , length );

  {
    typename vector_type::HostMirror hM =
      KokkosArray::create_mirror( matrix.values );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    KokkosArray::deep_copy( matrix.values , hM );
  }

  //------------------------------

  typename vector_type::HostMirror hx = KokkosArray::create_mirror( x );

  for ( size_t i = 0 ; i < length ; ++i ) {
    hx(i) = 1 + i ;
  }

  KokkosArray::deep_copy( x , hx );

  //------------------------------

  KokkosArray::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    KokkosArray::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename vector_type::HostMirror hy = KokkosArray::create_mirror( y );

    KokkosArray::deep_copy( hy , y );

    std::cout << std::endl << "test_flat_matrix"
              << std::endl ;
    for ( size_t i = 0 ; i < length ; ++i ) {
      std::cout << "hy(," << i << ") = " << hy(i) << std::endl ;
    }
  }

  return std::pair<size_t,double>( length , seconds_per_iter );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class Device >
void performance_test_driver_all( const int pdeg ,
				  const int minvar ,
				  const int maxvar ,
				  const int nGrid ,
				  const int nIter ,
				  const bool print ,
				  const bool test_block )
{
  typedef KokkosArray::NormalizedLegendrePolynomialBases<8> polynomial ;
  typedef KokkosArray::StochasticProductTensor< double , polynomial , Device , KokkosArray::CrsProductTensor > tensor_type ;

  std::cout.precision(8);

  //------------------------------

   std::vector< std::vector<size_t> > fem_graph ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );
  std::cout << std::endl << "\"FEM NNZ = " << graph_length << "\"" << std::endl;

  std::cout << std::endl
	    << "\"#nGrid\" , "
            << "\"#Variable\" , \"PolyDegree\" , \"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
	    << "\"Original-Flat MXV-Time\" , "
	    << "\"Original-Flat MXV-Speedup\" , "
            << "\"Original-Flat MXV-GFLOPS\" , "
	    << "\"Commuted-Flat MXV-Speedup\" , "
            << "\"Commuted-Flat MXV-GFLOPS\" , "
    //<< "\"Commuted-Flat NNZ\" , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
            << "\"Original-Matrix-Free-Block MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block MXV-GFLOPS\" , "
#endif
	    << "\"Block-Diagonal MXV-Speedup\" , "
            << "\"Block-Diagonal MXV-GFLOPS\" , "
    //<< "\"Block-Coord-Tensor MXV-Speedup\" , "
	    << "\"Block-Crs-Tensor MXV-Speedup\" , "
	    << "\"Block-Crs-Tensor MXV-Time\" , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
#endif
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const tensor_type tensor = KokkosArray::create_product_tensor< tensor_type >( var_degree );

    const std::vector<double> perf_matrix =
      test_product_tensor_diagonal_matrix<double,Device>( var_degree , nGrid , nIter , print );

    // const std::pair<size_t,double> perf_tensor =
    //   test_product_tensor_matrix<double,Device,KokkosArray::SparseProductTensor>( var_degree , nGrid , nIter , print );

    const std::pair<size_t,double> perf_crs_tensor =
      test_product_tensor_matrix<double,Device,KokkosArray::CrsProductTensor>( var_degree , nGrid , nIter , print );

    const std::vector<double> perf_flat_commuted =
      test_product_flat_commuted_matrix<double,Device>( var_degree , nGrid , nIter , print );

    const std::vector<double> perf_flat_original =
      test_product_flat_original_matrix<double,Device>( var_degree , nGrid , nIter , print );

#ifdef HAVE_KOKKOSARRAY_STOKHOS
    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<double,Device>( var_degree , nGrid , nIter , print , test_block );
#endif

    std::cout << nGrid << " , " << nvar << " , " << pdeg << " , "
	      << tensor.dimension() << " , "
	      << tensor.tensor().entry_count() << " , "
	      << perf_flat_original[0] << " , "
	      << perf_flat_original[1] << " , "
	      << perf_flat_original[1] / perf_flat_original[1] << " , "
              << perf_flat_original[2] << " , "
	      << perf_flat_original[1] / perf_flat_commuted[1] << " , "
              << perf_flat_commuted[2] << " , "
      //<< perf_flat_commuted[3] << " , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
	      << perf_flat_original[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
#endif
	      << perf_flat_original[1] / perf_matrix[1] << " , "
              << perf_matrix[2] << " , "
      //<< perf_flat_original.second / perf_tensor.second << " , "
	      << perf_flat_original[1] / perf_crs_tensor.second << " , "
	      << perf_crs_tensor.second << " , "
#ifdef HAVE_KOKKOSARRAY_STOKHOS
	      << perf_original_mat_free_block[3] / perf_crs_tensor.second << " , "
#endif
	      << std::endl ;
  }

  //------------------------------
}

#ifdef HAVE_KOKKOSARRAY_STOKHOS
template< class Device >
void performance_test_driver_poly( const int pdeg ,
				   const int minvar ,
				   const int maxvar ,
				   const int nGrid ,
				   const int nIter ,
				   const bool print ,
				   const bool test_block)
{
  typedef KokkosArray::NormalizedLegendrePolynomialBases<8> polynomial ;
  typedef KokkosArray::StochasticProductTensor< double , polynomial , Device , KokkosArray::CrsProductTensor > tensor_type ;

  std::cout.precision(8);

  //------------------------------

  std::vector< std::vector<size_t> > fem_graph ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );
  std::cout << std::endl << "\"FEM NNZ = " << graph_length << "\"" << std::endl;

  std::cout << std::endl
	    << "\"#nGrid\" , "
            << "\"#Variable\" , \"PolyDegree\" , \"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
	    << "\"Original-Matrix-Free-Block-MXV-Time\" , "
	    << "\"Original-Matrix-Free-Block-MXV-Speedup\" , "
            << "\"Original-Matrix-Free-Block-MXV-GFLOPS\" , "
	    << "\"Block-Crs-Tensor MXV-Time\" , "
	    << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-GFLOPS\" , "
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const tensor_type tensor = KokkosArray::create_product_tensor< tensor_type >( var_degree );

    const std::pair<size_t,double> perf_crs_tensor =
      test_product_tensor_matrix<double,Device,KokkosArray::CrsProductTensor>( var_degree , nGrid , nIter , print );

    const std::vector<double> perf_original_mat_free_block =
      test_original_matrix_free_vec<double,Device>( var_degree , nGrid , nIter , print , test_block );

    std::cout << nGrid << " , "
	      << nvar << " , " << pdeg << " , "
	      << tensor.dimension() << " , "
	      << tensor.tensor().entry_count() << " , "
	      << perf_original_mat_free_block[0] << " , "
	      << perf_original_mat_free_block[1] << " , "
	      << perf_original_mat_free_block[1] / perf_original_mat_free_block[1] << " , "
              << perf_original_mat_free_block[2] << " , "
	      << perf_crs_tensor.second << " , "
	      << perf_original_mat_free_block[1] / perf_crs_tensor.second << " , "
              << perf_original_mat_free_block[3] / perf_crs_tensor.second << " , "

	      << std::endl ;
  }

  //------------------------------
}
#endif

template< class Device >
void performance_test_driver(bool test_flat, bool test_orig, bool test_block)
{
}

//----------------------------------------------------------------------------

}


