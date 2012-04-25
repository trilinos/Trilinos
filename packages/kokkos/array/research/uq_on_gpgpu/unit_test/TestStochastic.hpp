
#include <utility>
#include <cmath>
#include <iostream>

namespace unit_test {

template< unsigned P >
void test_integration()
{
  Kokkos::GaussLegendre<P> rule ;
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

  Kokkos::GaussLegendre<P*2> rule ;
  Kokkos::NormalizedLegendrePolynomialBases<P,Device> poly ;

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

  Kokkos::GaussLegendre<P*3+1> rule ;
  Kokkos::NormalizedLegendrePolynomialBases<P,Device> poly ;

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
  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Device> polynomial ;
  typedef Kokkos::StochasticProductTensor< double , polynomial , Device , TensorType > tensor_type ;

  tensor_type tensor = Kokkos::create_product_tensor< tensor_type >( var_degree );

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
  typedef Kokkos::MultiVector< value_type , Device > vector_type ;

  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Device> polynomial ;

  typedef Kokkos::StochasticProductTensor< value_type , polynomial , Device , TensorType > tensor_type ;

  typedef Kokkos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;
  typedef typename matrix_type::graph_type graph_type ;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  const size_t graph_length = unit_test::generate_fem_graph( nGrid , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  matrix.block = Kokkos::create_product_tensor< tensor_type >( var_degree );
  matrix.graph = Kokkos::create_crsarray<graph_type>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();
  const size_t inner_matrix_size = matrix.block.dimension();

  matrix.values = Kokkos::create_multivector<vector_type>( inner_matrix_size , graph_length );

  vector_type x = Kokkos::create_multivector<vector_type>( inner_length , outer_length );
  vector_type y = Kokkos::create_multivector<vector_type>( inner_length , outer_length );

  typename vector_type::HostMirror hM = Kokkos::create_mirror( matrix.values );
  
  for ( size_t i = 0 ; i < graph_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1 + i ;
    }
  }
  
  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------
  // Generate input multivector:
  
  typename Kokkos::MultiVector<value_type,Device>::HostMirror hx = Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Kokkos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename Kokkos::MultiVector<value_type,Device>::HostMirror hy = Kokkos::create_mirror( y );

    Kokkos::deep_copy( hy , y );

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
std::pair<size_t,double>
test_product_tensor_diagonal_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef Kokkos::MultiVector<value_type,Device> vector_type ;

  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Kokkos::Host> polynomial ;
  typedef Kokkos::StochasticProductTensor< value_type , polynomial , Kokkos::Host , Kokkos::SparseProductTensor > tensor_type ;

  //------------------------------

  typedef Kokkos::Impl::Multiply<
            typename tensor_type::tensor_type ,
            Kokkos::SymmetricDiagonalSpec< Kokkos::Host > ,
            void > multiply_type ;

  typedef Kokkos::BlockCrsMatrix< Kokkos::SymmetricDiagonalSpec< Device > ,
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
    Kokkos::create_product_tensor< tensor_type >( var_degree );

  const size_t inner_length = tensor.dimension();

  //------------------------------
  // Generate CRS matrix of blocks with symmetric diagonal storage

  matrix_type matrix ;

  matrix.block  = Kokkos::SymmetricDiagonalSpec< Device >( inner_length );
  matrix.graph  = Kokkos::create_crsarray<graph_type>( std::string("test product tensor graph") , graph );
  matrix.values = Kokkos::create_multivector<vector_type>( matrix.block.matrix_size() , graph_length );

  Kokkos::MultiVector<value_type,Device> x = Kokkos::create_multivector<vector_type>( inner_length , outer_length );
  Kokkos::MultiVector<value_type,Device> y = Kokkos::create_multivector<vector_type>( inner_length , outer_length );

  typename Kokkos::MultiVector< value_type , Device >::HostMirror hM = Kokkos::create_mirror( matrix.values );

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

  Kokkos::deep_copy( matrix.values , hM );

  //------------------------------

  typename Kokkos::MultiVector< value_type , Device >::HostMirror hx = Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hx(j,i) = 1 + j + 10 * i ;
    }
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Kokkos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename Kokkos::MultiVector< value_type , Device >::HostMirror hy = Kokkos::create_mirror( y );

    Kokkos::deep_copy( hy , y );

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

  return std::pair<size_t,double>( outer_length * inner_length , seconds_per_iter );
}

//----------------------------------------------------------------------------
// Flatten to a plain CRS matrix

template< typename ScalarType , class Device >
std::pair<size_t,double>
test_product_flat_commuted_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef Kokkos::MultiVector< value_type , Device > vector_type ;

  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Kokkos::Host> polynomial ;

  typedef Kokkos::StochasticProductTensor<
     value_type , polynomial ,
     Kokkos::Host ,
     Kokkos::CrsProductTensor > tensor_type ;

  //------------------------------

  typedef Kokkos::CrsMatrix<value_type,Device> matrix_type ;
  typedef Kokkos::CrsArray<int,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const tensor_type tensor =
    Kokkos::create_product_tensor< tensor_type >( var_degree );

  const size_t inner_length = tensor.dimension();

  std::vector< std::vector<size_t> > tensor_graph( inner_length );

  for ( size_t iInnerRow = 0 ; iInnerRow < inner_length ; ++iInnerRow ) {
  for ( size_t iInnerCol = 0 ; iInnerCol < inner_length ; ++iInnerCol ) {

    for ( Kokkos::Host::size_type
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

  matrix.graph = Kokkos::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values =
    Kokkos::create_multivector<vector_type>( flat_graph_length );

  Kokkos::MultiVector<value_type,Device> x =
    Kokkos::create_multivector<vector_type>( flat_length );

  Kokkos::MultiVector<value_type,Device> y =
    Kokkos::create_multivector<vector_type>( flat_length );

  {
    typename Kokkos::MultiVector< value_type , Device >::HostMirror hM =
      Kokkos::create_mirror( matrix.values );

    for ( size_t i = 0 ; i < flat_graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    Kokkos::deep_copy( matrix.values , hM );
  }

  //------------------------------

  typename Kokkos::MultiVector< value_type , Device >::HostMirror hx =
    Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < flat_length ; ++i ) {
    hx(i) = 1 + i ;
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Kokkos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename Kokkos::MultiVector< value_type , Device >::HostMirror hy = Kokkos::create_mirror( y );

    Kokkos::deep_copy( hy , y );

    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
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
// Flatten to a plain CRS matrix

template< typename ScalarType , class Device >
std::pair<size_t,double>
test_product_flat_original_matrix(
  const std::vector<int> & var_degree ,
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef Kokkos::MultiVector< value_type , Device > vector_type ;

  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Kokkos::Host> polynomial ;

  typedef Kokkos::StochasticProductTensor<
    value_type , polynomial ,
    Kokkos::Host ,
    Kokkos::CrsProductTensor > tensor_type ;

  //------------------------------

  typedef Kokkos::CrsMatrix<value_type,Device> matrix_type ;
  typedef Kokkos::CrsArray<int,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t outer_length = nGrid * nGrid * nGrid ;
  unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const tensor_type tensor =
    Kokkos::create_product_tensor< tensor_type >( var_degree );

  const size_t inner_length = tensor.dimension();

  std::vector< std::vector<size_t> > tensor_graph( inner_length );

  for ( size_t iInnerRow = 0 ; iInnerRow < inner_length ; ++iInnerRow ) {
  for ( size_t iInnerCol = 0 ; iInnerCol < inner_length ; ++iInnerCol ) {

    for ( Kokkos::Host::size_type
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

  matrix.graph = Kokkos::create_crsarray<crsarray_type>( std::string("testing") , flat_graph );

  const size_t flat_graph_length = matrix.graph.entries.dimension(0);

  matrix.values =
    Kokkos::create_multivector<vector_type>( flat_graph_length );

  vector_type x = Kokkos::create_multivector<vector_type>( flat_length );
  vector_type y = Kokkos::create_multivector<vector_type>( flat_length );

  {
    typename Kokkos::MultiVector< value_type , Device >::HostMirror hM =
      Kokkos::create_mirror( matrix.values );

    for ( size_t i = 0 ; i < flat_graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    Kokkos::deep_copy( matrix.values , hM );
  }

  //------------------------------

  typename Kokkos::MultiVector< value_type , Device >::HostMirror hx =
    Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < flat_length ; ++i ) {
    hx(i) = 1 + i ;
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Kokkos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename Kokkos::MultiVector< value_type , Device >::HostMirror hy = Kokkos::create_mirror( y );

    Kokkos::deep_copy( hy , y );

    std::cout << std::endl << "test_product_flat_matrix"
              << std::endl ;
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
// A plain CRS matrix

template< typename ScalarType , class Device >
std::pair<size_t,double>
test_flat_matrix(
  const int nGrid ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;
  typedef Kokkos::MultiVector< value_type , Device > vector_type ;

  //------------------------------

  typedef Kokkos::CrsMatrix<value_type,Device> matrix_type ;
  typedef Kokkos::CrsArray<int,Device,int> crsarray_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > fem_graph ;

  const size_t length = nGrid * nGrid * nGrid ;
  const size_t graph_length =
    unit_test::generate_fem_graph( nGrid , fem_graph );

  //------------------------------

  matrix_type matrix ;

  matrix.graph = Kokkos::create_crsarray<crsarray_type>( std::string("testing") , fem_graph );

  matrix.values = Kokkos::create_multivector<vector_type>( graph_length );
  vector_type x = Kokkos::create_multivector<vector_type>( length );
  vector_type y = Kokkos::create_multivector<vector_type>( length );

  {
    typename Kokkos::MultiVector< value_type , Device >::HostMirror hM =
      Kokkos::create_mirror( matrix.values );

    for ( size_t i = 0 ; i < graph_length ; ++i ) {
      hM(i) = 1 + i ;
    }

    Kokkos::deep_copy( matrix.values , hM );
  }

  //------------------------------

  typename Kokkos::MultiVector< value_type , Device >::HostMirror hx =
    Kokkos::create_mirror( x );

  for ( size_t i = 0 ; i < length ; ++i ) {
    hx(i) = 1 + i ;
  }

  Kokkos::deep_copy( x , hx );

  //------------------------------

  Kokkos::Impl::Timer clock ;
  for ( int iter = 0 ; iter < iterCount ; ++iter ) {
    Kokkos::multiply( matrix , x , y );
  }
  Device::fence();

  const double seconds_per_iter = clock.seconds() / ((double) iterCount );

  //------------------------------

  if ( print_flag ) {
    typename Kokkos::MultiVector< value_type , Device >::HostMirror hy = Kokkos::create_mirror( y );

    Kokkos::deep_copy( hy , y );

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
void performance_test_driver( const int pdeg ,
                              const int minvar ,
                              const int maxvar ,
                              const int nGrid ,
                              const int nIter ,
                              const bool print )
{
  typedef Kokkos::NormalizedLegendrePolynomialBases<8,Device> polynomial ;
  typedef Kokkos::StochasticProductTensor< double , polynomial , Device , Kokkos::SparseProductTensor > tensor_type ;

  std::cout.precision(8);

  //------------------------------

  std::cout << std::endl
            << "\"#Variable\" , \"PolyDegree\" , \"#Bases\" , "
            << "\"#TensorEntry\" , "
            << "\"VectorSize\" , "
            << "\"Original-Flat MXV-Time\" , "
            << "\"Original-Flat MXV-Speedup\" , "
            << "\"Commuted-Flat MXV-Speedup\" , "
            << "\"Block-Coord-Tensor MXV-Speedup\" , "
            << "\"Block-Crs-Tensor MXV-Speedup\" , "
            << "\"Block-Diagonal MXV-Speedup\""
            << std::endl ;

  for ( int nvar = minvar ; nvar <= maxvar ; ++nvar ) {
    std::vector<int> var_degree( nvar , pdeg );

    const tensor_type tensor = Kokkos::create_product_tensor< tensor_type >( var_degree );

    const std::pair<size_t,double> perf_matrix =
      test_product_tensor_diagonal_matrix<double,Device>( var_degree , nGrid , nIter , print );

    const std::pair<size_t,double> perf_tensor =
      test_product_tensor_matrix<double,Device,Kokkos::SparseProductTensor>( var_degree , nGrid , nIter , print );

    const std::pair<size_t,double> perf_crs_tensor =
      test_product_tensor_matrix<double,Device,Kokkos::CrsProductTensor>( var_degree , nGrid , nIter , print );

    const std::pair<size_t,double> perf_flat_commuted =
      test_product_flat_commuted_matrix<double,Device>( var_degree , nGrid , nIter , print );

    const std::pair<size_t,double> perf_flat_original =
      test_product_flat_original_matrix<double,Device>( var_degree , nGrid , nIter , print );

     std::cout << nvar << " , " << pdeg << " , "
               << tensor.dimension() << " , "
               << tensor.tensor().entry_count() << " , "
               << perf_flat_original.first << " , "
               << perf_flat_original.second << " , "
               << perf_flat_original.second / perf_flat_original.second << " , "
               << perf_flat_original.second / perf_flat_commuted.second << " , "
               << perf_flat_original.second / perf_tensor.second << " , "
               << perf_flat_original.second / perf_crs_tensor.second << " , "
               << perf_flat_original.second / perf_matrix.second
               << std::endl ;
  }

  //------------------------------
}

template< class Device >
void performance_test_driver()
{
  const int nGrid = 5 ;
  const int nIter = 10 ; 
  const bool print = false ;

  performance_test_driver<Device>( 3 , 4 , 11 , nGrid , nIter , print );
  performance_test_driver<Device>( 4 , 3 ,  7 , nGrid , nIter , print );

  //------------------------------

  std::cout << std::endl
            << "\"CRS flat-matrix ~27 nonzeros/row (CUDA uses cusparse)\""
            << std::endl
            << "\"VectorSize\" , "
            << "\"MXV-Time\""
            << std::endl ;

  for ( int n_grid = nGrid ; n_grid <= ( nGrid << 5 ) ; n_grid <<= 1 ) {

    const std::pair<size_t,double> perf_flat =
      test_flat_matrix<double,Device>( n_grid , nIter , print );

    std::cout << perf_flat.first << " , "
              << perf_flat.second
              << std::endl ;
  }

  //------------------------------
}

//----------------------------------------------------------------------------

}


