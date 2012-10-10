
#include <KokkosEmbed_Array.hpp>
#include <impl/KokkosEmbed_Array_ViewOperRight.hpp>
#include <iostream>

template< class Device > int test();

template<>
int test< TEST_KOKKOSARRAY_DEVICE >()
{
  KokkosArray::Array<double,10> a = 1 ;
  KokkosArray::Array<double,10> b = a + a ;

  volatile KokkosArray::Array<double,10> c = a ;

  for ( unsigned i = 0 ; i < 10 ; ++i ) a[i] = i + 1 ;

  std::cout << -a << std::endl ;
  std::cout << a + b << std::endl ;
  std::cout << b - a << std::endl ;
  std::cout << b / a << std::endl ;
  std::cout << b / ( b + a ) << std::endl ;

  b = 7 ;
  b *= b ;
  std::cout << sqrt( b ) << std::endl ;
  std::cout << cos( 1 + sqrt( a ) ) << std::endl ;

  std::cout << 5 * ( b + a ) / ( b - a ) << std::endl ;

  std::cout << ( b = c * ( a + c ) / 10 ) << std::endl ;

  KokkosArray::print_expression( std::cout , a );
  KokkosArray::print_expression( std::cout , 5 * ( b + a ) / ( b - a ) );
  KokkosArray::print_expression( std::cout , a * c + b );

  // b + b = c ; // test compilation error

  int x_raw[10] ;
  int y_raw[100] ;

  KokkosArray::Array<int,10,KokkosArray::ArrayProxyContiguous> X( x_raw );
  KokkosArray::Array<int,10,KokkosArray::ArrayProxyStrided> Y( y_raw , 10 );

  Y = b ;
  X = a ;
  X = abs( X * Y );
  X = max( X , Y );
  X = min( X , Y * 2 );
  X = min( X , 2 );
  X = min( 2 , X );

  KokkosArray::print_expression( std::cout , min( X , Y * 2 ) );

  a = X ;

  X = 2 * Y ;

  if ( X < Y ) {
    std::cout << X << " < " << Y << std::endl ;
  }
  else if ( Y < X ) {
    std::cout << Y << " < " << X << std::endl ;
  }
  else if ( X == Y ) {
    std::cout << Y << " == " << X << std::endl ;
  }
  else {
    std::cout << Y << " NotOrdered " << X << std::endl ;
  }

  if ( b < X ) {
    std::cout << b << " < " << X << std::endl ;
  }
  else if ( b > X ) {
    std::cout << b << " > " << X << std::endl ;
  }
  else if ( b == X ) {
    std::cout << b << " == " << X << std::endl ;
  }
  else {
    std::cout << b << " NotOrdered " << X << std::endl ;
  }

  Y = X * b ;

  return 0 ;
}

