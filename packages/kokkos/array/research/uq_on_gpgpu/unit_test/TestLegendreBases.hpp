
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
  typedef Kokkos::StochasticProductTensor< double , polynomial , Device > bases_type ;

  bases_type bases = Kokkos::create_product_tensor< bases_type >( var_degree );

  // Verification?
}

//----------------------------------------------------------------------------

template< class ScalarType , class Device >
void test_product_tensor_matrix(
  const std::vector<int> & var_degree ,
  const int nGraph )
{
  typedef ScalarType value_type ;

  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Device> polynomial ;

  typedef Kokkos::StochasticProductTensor< value_type , polynomial , Device > bases_type ;

  typedef Kokkos::BlockCrsMatrix< bases_type , value_type , Device > matrix_type ;

  std::vector< std::vector<size_t> > graph ;

  const size_t outer_length = unit_test::generate_fem_graph( nGraph , graph );

  matrix_type matrix ;

  matrix.block = Kokkos::create_product_tensor< bases_type >( var_degree );
  matrix.graph = Kokkos::create_labeled_crsmap<Device>( std::string("test crs graph") , graph );

  const size_t inner_length = matrix.block.dimension();
  const size_t inner_matrix_size = matrix.block.dimension();

  matrix.values = Kokkos::create_multivector<value_type,Device>( inner_matrix_size , outer_length );

  Kokkos::MultiVector<value_type,Device> x = Kokkos::create_multivector<value_type,Device>( inner_length , outer_length );
  Kokkos::MultiVector<value_type,Device> y = Kokkos::create_multivector<value_type,Device>( inner_length , outer_length );

  typename Kokkos::MultiVector<value_type,Device>::HostMirror hM = Kokkos::create_mirror( matrix.values );
  typename Kokkos::MultiVector<value_type,Device>::HostMirror hx = Kokkos::create_mirror( x );
  typename Kokkos::MultiVector<value_type,Device>::HostMirror hy = Kokkos::create_mirror( y );

  for ( size_t i = 0 ; i < outer_length ; ++i ) {
    for ( size_t j = 0 ; j < inner_length ; ++j ) {
      hM(j,i) = 1 + i ;
      hx(j,i) = 1 + j + 10 * i ;
    }
  }
  
  Kokkos::deep_copy( x , hx );
  Kokkos::deep_copy( matrix.values , hM );
  
  Kokkos::multiply( matrix , x , y );
  
  Kokkos::deep_copy( hy , y );
}


}


