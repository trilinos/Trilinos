
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

template< class Device >
void test_product_tensor( const std::vector<int> & var_degree )
{
  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Device> polynomial ;
  typedef Kokkos::StochasticProductTensor< double , polynomial , Device > tensor_type ;

  tensor_type tensor = Kokkos::create_product_tensor< tensor_type >( var_degree );

  // Verification?
}

//----------------------------------------------------------------------------

template< typename ScalarType , class Device >
std::pair<double,double>
test_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGraph ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;

  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Device> polynomial ;

  typedef Kokkos::StochasticProductTensor< value_type , polynomial , Device > tensor_type ;

  typedef Kokkos::BlockCrsMatrix< tensor_type , value_type , Device > matrix_type ;

  //------------------------------
  // Generate graph for "FEM" box structure:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGraph * nGraph * nGraph ;
  const size_t graph_length = unit_test::generate_fem_graph( nGraph , graph );

  //------------------------------
  // Generate CRS block-tensor matrix:

  matrix_type matrix ;

  matrix.block = Kokkos::create_product_tensor< tensor_type >( var_degree );
  matrix.graph = Kokkos::create_labeled_crsmap<Device>( std::string("test crs graph") , graph );

  const size_t inner_length      = matrix.block.dimension();
  const size_t inner_matrix_size = matrix.block.dimension();

  matrix.values = Kokkos::create_multivector<value_type,Device>( inner_matrix_size , graph_length );

  Kokkos::MultiVector<value_type,Device> x = Kokkos::create_multivector<value_type,Device>( inner_length , outer_length );
  Kokkos::MultiVector<value_type,Device> y = Kokkos::create_multivector<value_type,Device>( inner_length , outer_length );

  typename Kokkos::MultiVector<value_type,Device>::HostMirror hM = Kokkos::create_mirror( matrix.values );
  
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
  const double seconds = clock.seconds();

  const double effective_flop_rate =
    ((size_t) iterCount ) *
    ((size_t) matrix.graph.entry_count() ) *
    ((size_t) matrix.block.dimension() ) *
    ((size_t) matrix.block.dimension() ) / seconds ;

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

  return std::pair<double,double>( seconds , effective_flop_rate );
}

//----------------------------------------------------------------------------

template< typename ScalarType , class Device >
std::pair<double,double>
test_product_tensor_diagonal_matrix(
  const std::vector<int> & var_degree ,
  const int nGraph ,
  const int iterCount ,
  const bool print_flag = false )
{
  typedef ScalarType value_type ;

  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Kokkos::Host> polynomial ;
  typedef Kokkos::StochasticProductTensor< value_type , polynomial , Kokkos::Host > tensor_type ;

  //------------------------------

  typedef Kokkos::Impl::Multiply<
            typename tensor_type::tensor_type ,
            Kokkos::SymmetricDiagonalSpec< Kokkos::Host > ,
            void > multiply_type ;

  typedef Kokkos::BlockCrsMatrix< Kokkos::SymmetricDiagonalSpec< Device > ,
                                  value_type , Device > matrix_type ;

  //------------------------------
  // Generate FEM graph:

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = nGraph * nGraph * nGraph ;
  const size_t graph_length = unit_test::generate_fem_graph( nGraph , graph );

  //------------------------------
  // Generate product tensor from variables' degrees

  const tensor_type tensor =
    Kokkos::create_product_tensor< tensor_type >( var_degree );

  const size_t inner_length = tensor.dimension();

  //------------------------------
  // Generate CRS matrix of blocks with symmetric diagonal storage

  matrix_type matrix ;

  matrix.block  = Kokkos::SymmetricDiagonalSpec< Device >( inner_length );
  matrix.graph  = Kokkos::create_labeled_crsmap<Device>( std::string("test product tensor graph") , graph );
  matrix.values = Kokkos::create_multivector<value_type,Device>( matrix.block.matrix_size() , graph_length );

  Kokkos::MultiVector<value_type,Device> x = Kokkos::create_multivector<value_type,Device>( inner_length , outer_length );
  Kokkos::MultiVector<value_type,Device> y = Kokkos::create_multivector<value_type,Device>( inner_length , outer_length );

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
  const double seconds = clock.seconds();

  const double effective_flop_rate =
    ((size_t) iterCount ) *
    ((size_t) matrix.graph.entry_count() ) *
    ((size_t) matrix.block.dimension() ) *
    ((size_t) matrix.block.dimension() ) / seconds ;

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

  return std::pair<double,double>( seconds , effective_flop_rate );
}

//----------------------------------------------------------------------------

}


