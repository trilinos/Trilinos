/*
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
*/

#include <stdio.h>
#include <iostream>
#include <typeinfo>

#include <Kokkos_Core.hpp>
#include <Kokkos_Array.hpp>

#include <impl/Kokkos_ArrayAnalyzeShape.hpp>
#include <impl/Kokkos_ArrayViewDefault.hpp>

//----------------------------------------------------------------------------

template< class Space > struct TestDevice ;

template<> struct TestDevice< Kokkos::HostSpace >
{
  enum { value = true };

  typedef Kokkos::Threads type ;

  TestDevice()
  {
    Kokkos::Threads::initialize( 4 , 1 /* use one NUMA region */ );
  }
  ~TestDevice()
  {
    Kokkos::Threads::finalize();
  }
};

template<>
struct TestDevice< Kokkos::CudaSpace >
{
#if defined( __CUDACC__ )

  enum { value = true };

  typedef Kokkos::Cuda type ;

  TestDevice()
  {
    Kokkos::Cuda::initialize();
  }
  ~TestDevice()
  {
    Kokkos::Cuda::finalize();
  }
#else
  enum { value = false };

  typedef Kokkos::Threads type ;
#endif
};

//----------------------------------------------------------------------------

template< class Device > int test();

template<>
int test< TEST_KOKKOS_SPACE >()
{
  if ( ! TestDevice< TEST_KOKKOS_SPACE >::value ) return 0 ;

  Kokkos::Array<double,10> a = 1 ;
  Kokkos::Array<double,10> b = a + a ;

  volatile Kokkos::Array<double,10> c = a ;

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

  Kokkos::print_expression( std::cout , a );
  Kokkos::print_expression( std::cout , 5 * ( b + a ) / ( b - a ) );
  Kokkos::print_expression( std::cout , a * c + b );

  // b + b = c ; // test compilation error

  int x_raw[10] ;
  int y_raw[100] ;

  Kokkos::Array<int,10,Kokkos::ArrayProxyContiguous> X( x_raw );
  Kokkos::Array<int,10,Kokkos::ArrayProxyStrided> Y( y_raw , 10 );

  Y = b ;
  X = a ;
  X = abs( X * Y );
  X = max( X , Y );
  X = min( X , Y * 2 );
  X = min( X , 2 );
  X = min( 2 , X );

  Kokkos::print_expression( std::cout , min( X , Y * 2 ) );

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

template< class Device > int testdyn();

template<>
int testdyn< TEST_KOKKOS_SPACE >()
{
  if ( ! TestDevice< TEST_KOKKOS_SPACE >::value ) return 0 ;

  enum { Count = 10 };

  Kokkos::Array<double,0> a(1, Count);
  Kokkos::Array<double,0> b = a + a ;

  volatile Kokkos::Array<double,0> c = a ;

  for ( unsigned i = 0 ; i < 10 ; ++i ) a[i] = i + 1 ;

  std::cout << a << std::endl ;
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

  Kokkos::print_expression( std::cout , a );
  Kokkos::print_expression( std::cout , 5 * ( b + a ) / ( b - a ) );
  Kokkos::print_expression( std::cout , a * c + b );

  // b + b = c ; // test compilation error

  int x_raw[10] ;
  int y_raw[100] ;

  Kokkos::Array<int,10,Kokkos::ArrayProxyContiguous> X( x_raw );
  Kokkos::Array<int,10,Kokkos::ArrayProxyStrided> Y( y_raw , 10 );

  Y = b ;
  X = a ;
  X = abs( X * Y );
  X = max( X , Y );
  X = min( X , Y * 2 );
  X = min( X , 2 );
  X = min( 2 , X );

  Kokkos::print_expression( std::cout , min( X , Y * 2 ) );

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
//----------------------------------------------------------------------------

template< class ValueType , class DeviceType >
struct TestFunctor
{
  typedef DeviceType                       execution_space ;
  typedef typename execution_space::size_type  size_type ;

  typedef Kokkos::View< ValueType * , execution_space > vector_type ;

  const vector_type x , y , z ;
  const ValueType   a , b ;

  explicit
  TestFunctor( const unsigned N )
    : x("x",N), y("y",N), z("z",N)
    , a(2), b(3)
    { }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type ip ) const
  {
    z(ip) = a * x(ip) + b * y(ip);
  }
};

//----------------------------------------------------------------------------

template< class Space > int test_functor();

template<>
int test_functor< TEST_KOKKOS_SPACE >()
{
  typedef TestDevice< TEST_KOKKOS_SPACE > TestD ;
  typedef TestD::type execution_space ;

  if ( ! TestD::value ) return 0 ;
  
  enum { Count = 20 };

  typedef Kokkos::Array<double,Count> value_type ;

  typedef TestFunctor<value_type,execution_space> functor_type ;

  TestD device ;

  std::cout << "functor_type::vector_type::scalar_type = "
            << typeid(functor_type::vector_type::scalar_type).name()
            << std::endl ;

  std::cout << "functor_type::vector_type::value_type = "
            << typeid(functor_type::vector_type::value_type).name()
            << std::endl ;

  const unsigned N = 1000 ;

  functor_type f( N );

  functor_type::vector_type::HostMirror x = 
    Kokkos::create_mirror( f.x );

  std::cout << "  distance( x(0)[0] , x(0)[last] ) = "
            << (int)( & x(0)[Count-1] - & x(0)[0] )
            << std::endl
            << "  x.shape() ="
            << " { " << x.shape().N0
            << " , " << x.shape().N1
            << " }" <<" }" <<  std::endl ;

  if ( & x(0)[0] != & x(0,0) ) {
    std::cout << "  FAILED & x(0)[0] != & x(0,0) : "
              << & x(0)[0] 
              << " != "
              << & x(0,0)
              << std::endl ;
  }

  for ( unsigned i = 0 ; i < N ; ++i ) { x(i) = 1 ; }

  Kokkos::deep_copy( f.x , x );
  Kokkos::deep_copy( f.y , x );

  Kokkos::parallel_for( N , f );

  Kokkos::deep_copy( x , f.z );

  std::cout << "test_functor z = " << x(0) << std::endl ;

  return 0 ;
}

//----------------------------------------------------------------------------

template< typename ScalarType , unsigned N , class ProxyX , class ProxyY >
KOKKOS_INLINE_FUNCTION
double dot( const Kokkos::Array< ScalarType , N , ProxyX > & x ,
            const Kokkos::Array< ScalarType , N , ProxyY > & y )
{
  double r = 0 ;
  for ( unsigned i = 0 ; i < x.size() ; ++i ) {
    r += x[i] * y[i] ;
  }
  return r ;
}

template< class ScalarType , class DeviceType >
struct TestDot
{
  typedef DeviceType                       execution_space ;
  typedef typename execution_space::size_type  size_type ;

  typedef double value_type ;

  typedef Kokkos::View< ScalarType * , execution_space > vector_type ;

  const vector_type x , y ;

  TestDot( const vector_type & arg_x ,
           const vector_type & arg_y )
    : x( arg_x ), y( arg_y ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type ip , value_type & update ) const
    { update += dot( x(ip) , y(ip) ); }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update += input ; }
};


template< class Space > int test_inner_product();

template<>
int test_inner_product< TEST_KOKKOS_SPACE >()
{
  typedef TestDevice< TEST_KOKKOS_SPACE > TestD ;
  typedef TestD::type execution_space ;

  if ( ! TestD::value ) return 0 ;

  TestD device ;

  enum { Count = 20 };

  typedef Kokkos::Array<double,Count> value_type ;

  typedef TestDot< value_type , execution_space > functor_type ;

  const unsigned N = 1000 ;

  functor_type::vector_type x("x",N);
  functor_type::vector_type y("y",N);

  functor_type::vector_type::HostMirror hx = Kokkos::create_mirror_view( x );
  functor_type::vector_type::HostMirror hy = Kokkos::create_mirror_view( y );
  
  for ( unsigned i = 0 ; i < N ; ++i ) { hx(i) = 1 ; hy(i) = 2 ; }

  Kokkos::deep_copy( x , hx );
  Kokkos::deep_copy( y , hy );

  double r = 0 ;

  Kokkos::parallel_reduce( N , functor_type(x,y) , r );

  std::cout << "test_inner_product r = " << r << std::endl ;

  return 0 ;
}


