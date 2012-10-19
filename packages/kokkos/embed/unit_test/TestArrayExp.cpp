#include <iostream>
#include <typeinfo>

#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_Host.hpp>
#include <KokkosArray_Array.hpp>

#include <impl/KokkosArray_ArrayAnalyzeShape.hpp>
#include <impl/KokkosArray_ArrayViewOperRight.hpp>
#include <impl/KokkosArray_ArrayViewOperLeft.hpp>

//----------------------------------------------------------------------------

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

//----------------------------------------------------------------------------

template< class Device > struct TestFunctor ;

template<>
struct TestFunctor< TEST_KOKKOSARRAY_DEVICE >
{
  typedef KokkosArray::Array<double,20> scalar_type ;

  typedef TEST_KOKKOSARRAY_DEVICE device_type ;
  typedef device_type::size_type  size_type ;

  typedef KokkosArray::View< scalar_type * , device_type > vector_type ;

  static const size_type N = 1000 ;

  const vector_type x , y , z ;
  const scalar_type a , b ;

  TestFunctor()
    : x("x",N), y("y",N), z("z",N)
    , a(2), b(3)
    { }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void operator()( const size_type ip ) const
  {
    z(ip) = a * x(ip) + b * y(ip);
  }
};

template< class Device > int test_functor();

template<>
int test_functor< TEST_KOKKOSARRAY_DEVICE >()
{
  typedef TEST_KOKKOSARRAY_DEVICE device_type ;
  
  typedef TestFunctor<device_type> functor_type ;

  std::cout << "functor_type::vector_type::scalar_type = "
            << typeid(functor_type::vector_type::scalar_type).name()
            << std::endl ;

  std::cout << "functor_type::vector_type::value_type = "
            << typeid(functor_type::vector_type::value_type).name()
            << std::endl ;

  functor_type f ;

  functor_type::vector_type::HostMirror x = 
    KokkosArray::create_mirror( f.x );

  for ( unsigned i = 0 ; i < f.N ; ++i ) { x(i) = 1 ; }

  KokkosArray::deep_copy( f.x , x );
  KokkosArray::deep_copy( f.y , x );

  KokkosArray::parallel_for( f.N , f );

  KokkosArray::deep_copy( x , f.z );

  std::cout << "test_functor z = " << x(0) << std::endl ;

  return 0 ;
}


