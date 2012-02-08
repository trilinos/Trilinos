
#include <cmath>
#include <iostream>

namespace Test {

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


template< class Device >
void test_product_tensor( const std::vector<int> & var_degree )
{
  typedef Kokkos::NormalizedLegendrePolynomialBases<4,Device> polynomial ;
  typedef Kokkos::ProductTensorFromBases< double , polynomial , Device > bases_type ;

  bases_type bases = Kokkos::create_product_tensor< bases_type >( var_degree );

  // Verification?
}


}


