/*
// ************************************************************************
//
//                             KokkosArray
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

#include <KokkosArray_Cuda.hpp>
#include <KokkosArray_Host.hpp>
#include <KokkosArray_Array.hpp>

#include <impl/KokkosArray_ArrayAnalyzeShape.hpp>
#include <impl/KokkosArray_ArrayViewOperRight.hpp>
#include <impl/KokkosArray_ArrayViewOperLeft.hpp>

//----------------------------------------------------------------------------

template< class Space > struct TestDevice ;

template<> struct TestDevice< KokkosArray::HostSpace >
{
  enum { value = true };

  typedef KokkosArray::Host type ;

  TestDevice()
  {
    KokkosArray::Host::initialize( 1 , 4 );
  }
  ~TestDevice()
  {
    KokkosArray::Host::finalize();
  }
};

template<>
struct TestDevice< KokkosArray::CudaSpace >
{
#if defined( __CUDACC__ )

  enum { value = true };

  typedef KokkosArray::Cuda type ;

  TestDevice()
  {
    KokkosArray::Cuda::initialize();
  }
  ~TestDevice()
  {
    KokkosArray::Cuda::finalize();
  }
#else
  enum { value = false };

  typedef KokkosArray::Host type ;
#endif
};

//----------------------------------------------------------------------------

template< class Device > int test();

template<>
int test< TEST_KOKKOSARRAY_SPACE >()
{
  if ( ! TestDevice< TEST_KOKKOSARRAY_SPACE >::value ) return 0 ;

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
//----------------------------------------------------------------------------

template< class ValueType , class DeviceType >
struct TestFunctor
{
  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;

  typedef KokkosArray::View< ValueType * , device_type > vector_type ;

  const vector_type x , y , z ;
  const ValueType   a , b ;

  explicit
  TestFunctor( const unsigned N )
    : x("x",N), y("y",N), z("z",N)
    , a(2), b(3)
    { }

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type ip ) const
  {
    z(ip) = a * x(ip) + b * y(ip);
  }
};

//----------------------------------------------------------------------------

template< class Space > int test_functor();

template<>
int test_functor< TEST_KOKKOSARRAY_SPACE >()
{
  typedef TestDevice< TEST_KOKKOSARRAY_SPACE > TestD ;
  typedef TestD::type device_type ;

  if ( ! TestD::value ) return 0 ;
  
  enum { Count = 20 };

  typedef KokkosArray::Array<double,Count> value_type ;

  typedef TestFunctor<value_type,device_type> functor_type ;

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
    KokkosArray::create_mirror( f.x );

  std::cout << "  distance( x(0)[0] , x(0)[last] ) = "
            << (int)( & x(0)[Count-1] - & x(0)[0] )
            << std::endl
            << "  x.shape() = { " << x.shape().Stride
            << " : " << x.shape().N0 
            << " , " << x.shape().N1
            << std::endl ;

  if ( & x(0)[0] != & x(0,0) ) {
    std::cout << "  FAILED & x(0)[0] != & x(0,0) : "
              << & x(0)[0] 
              << " != "
              << & x(0,0)
              << std::endl ;
  }

  for ( unsigned i = 0 ; i < N ; ++i ) { x(i) = 1 ; }

  KokkosArray::deep_copy( f.x , x );
  KokkosArray::deep_copy( f.y , x );

  KokkosArray::parallel_for( N , f );

  KokkosArray::deep_copy( x , f.z );

  std::cout << "test_functor z = " << x(0) << std::endl ;

  return 0 ;
}

//----------------------------------------------------------------------------

template< typename ScalarType , unsigned N , class ProxyX , class ProxyY >
KOKKOSARRAY_INLINE_FUNCTION
double dot( const KokkosArray::Array< ScalarType , N , ProxyX > & x ,
            const KokkosArray::Array< ScalarType , N , ProxyY > & y )
{
  double r = 0 ;
  for ( unsigned i = 0 ; i < x.value_count ; ++i ) {
    r += x[i] * y[i] ;
  }
  return r ;
}

template< class ScalarType , class DeviceType >
struct TestDot
{
  typedef DeviceType                       device_type ;
  typedef typename device_type::size_type  size_type ;

  typedef double value_type ;

  typedef KokkosArray::View< ScalarType * , device_type > vector_type ;

  const vector_type x , y ;

  TestDot( const vector_type & arg_x ,
           const vector_type & arg_y )
    : x( arg_x ), y( arg_y ) {}

  KOKKOSARRAY_INLINE_FUNCTION
  void operator()( const size_type ip , value_type & update ) const
    { update += dot( x(ip) , y(ip) ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update += input ; }
};


template< class Space > int test_inner_product();

template<>
int test_inner_product< TEST_KOKKOSARRAY_SPACE >()
{
  typedef TestDevice< TEST_KOKKOSARRAY_SPACE > TestD ;
  typedef TestD::type device_type ;

  if ( ! TestD::value ) return 0 ;

  TestD device ;

  enum { Count = 20 };

  typedef KokkosArray::Array<double,Count> value_type ;

  typedef TestDot< value_type , device_type > functor_type ;

  const unsigned N = 1000 ;

  functor_type::vector_type x("x",N);
  functor_type::vector_type y("y",N);

  functor_type::vector_type::HostMirror hx = KokkosArray::create_mirror_view( x );
  functor_type::vector_type::HostMirror hy = KokkosArray::create_mirror_view( y );
  
  for ( unsigned i = 0 ; i < N ; ++i ) { hx(i) = 1 ; hy(i) = 2 ; }

  KokkosArray::deep_copy( x , hx );
  KokkosArray::deep_copy( y , hy );

  double r = KokkosArray::parallel_reduce( N , functor_type(x,y) );

  std::cout << "test_inner_product r = " << r << std::endl ;

  return 0 ;
}


