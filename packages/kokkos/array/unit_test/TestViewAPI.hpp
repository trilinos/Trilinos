/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <impl/KokkosArray_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace {

template< typename T, class > class TestViewAPI ;
template< typename T, class > class TestViewOperator ;

template< typename T >
struct TestViewOperator< T , KOKKOS_MACRO_DEVICE >
{
  typedef KOKKOS_MACRO_DEVICE  device_type ;

  static const unsigned N = 100 ;
  static const unsigned D = 3 ;


  typedef KokkosArray::View< T[][D] , device_type > view_type ;

  const view_type v1 ;
  const view_type v2 ;

  TestViewOperator()
    : v1( KokkosArray::create<view_type>( "v1" , N ) )
    , v2( KokkosArray::create<view_type>( "v1" , N ) )
    {}

  static void apply()
  {
    KokkosArray::parallel_for( N , TestViewOperator() );
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const unsigned i ) const
  {
    const unsigned X = 0 ;
    const unsigned Y = 1 ;
    const unsigned Z = 2 ;

    v2(i,X) = v1(i,X);
    v2(i,Y) = v1(i,Y);
    v2(i,Z) = v1(i,Z);
  }
};


template<typename T>
class TestViewAPI< T, KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE  device ;
  typedef KokkosArray::Host    host ;

  TestViewAPI()
  {
    run_test();
    run_test_const();
    run_test_vector();
    TestViewOperator< T , device >::apply();
  }

  enum { NP = 1000 ,
         N1 = 3 ,
         N2 = 5 ,
         N3 = 7 };

  typedef KokkosArray::View< T , device > dView0 ;
  typedef KokkosArray::View< T[] , device > dView1 ;
  typedef KokkosArray::View< T[][N1] , device > dView2 ;
  typedef KokkosArray::View< T[][N1][N2] , device > dView3 ;
  typedef KokkosArray::View< T[][N1][N2][N3] , device > dView4 ;
  typedef KokkosArray::View< const T[][N1][N2][N3] , device > const_dView4 ;

  struct MyView4 {
    typedef T       data_type[][N1][N2][N3] ;
    typedef device  layout_type ;
    typedef device  device_type ;
  };

  static void run_test()
  {
    typedef typename dView0::HostMirror  hView0 ;
    typedef typename dView1::HostMirror  hView1 ;
    typedef typename dView2::HostMirror  hView2 ;
    typedef typename dView3::HostMirror  hView3 ;
    typedef typename dView4::HostMirror  hView4 ;

    dView4 dx , dy , dz ;
    hView4 hx , hy , hz ;
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_FALSE(dz);
    ASSERT_FALSE(hx);
    ASSERT_FALSE(hy);
    ASSERT_FALSE(hz);

    dx = KokkosArray::create<dView4>( "dx" , NP );
    dy = KokkosArray::create<MyView4>( "dy" , NP );

    const_dView4 const_dx = dx ;

    ASSERT_TRUE(dx);
    ASSERT_TRUE(const_dx);
    ASSERT_TRUE(dy);
    ASSERT_NE( dx , dy );

    ASSERT_EQ( dx.dimension(0) , NP );
    ASSERT_EQ( dx.dimension(1) , N1 );
    ASSERT_EQ( dx.dimension(2) , N2 );
    ASSERT_EQ( dx.dimension(3) , N3 );

    ASSERT_EQ( dy.dimension(0) , NP );
    ASSERT_EQ( dy.dimension(1) , N1 );
    ASSERT_EQ( dy.dimension(2) , N2 );
    ASSERT_EQ( dy.dimension(3) , N3 );

    hx = KokkosArray::create_mirror( dx );
    hy = KokkosArray::create_mirror( dy );

    size_t count = 0 ;
    for ( size_t ip = 0 ; ip < NP ; ++ip ) {
    for ( size_t i1 = 0 ; i1 < hx.dimension(1) ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < hx.dimension(2) ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < hx.dimension(3) ; ++i3 ) {
      hx(ip,i1,i2,i3) = ++count ;
    }}}}

    KokkosArray::deep_copy( dx , hx );
    KokkosArray::deep_copy( dy , dx );
    KokkosArray::deep_copy( hy , dy );

    for ( size_t ip = 0 ; ip < NP ; ++ip ) {
    for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
      { ASSERT_EQ( hx(ip,i1,i2,i3) , hy(ip,i1,i2,i3) ); }
    }}}}

    dz = dx ; ASSERT_EQ( dx, dz); ASSERT_NE( dy, dz);
    dz = dy ; ASSERT_EQ( dy, dz); ASSERT_NE( dx, dz);

    dx = dView4();
    ASSERT_FALSE(dx);
    ASSERT_TRUE( dy);
    ASSERT_TRUE( dz);
    dy = dView4();
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_TRUE( dz);
    dz = dView4();
    ASSERT_FALSE(dx);
    ASSERT_FALSE(dy);
    ASSERT_FALSE(dz);
  }

  typedef T DataType[2] ;

  static void
  check_auto_conversion_to_const(
     const KokkosArray::View< const DataType , device > & arg_const ,
     const KokkosArray::View< DataType , device > & arg )
  {
    ASSERT_TRUE( arg_const == arg );
  }

  static void run_test_const()
  {
    typedef KokkosArray::View< DataType , device > typeX ;
    typedef KokkosArray::View< const DataType , device > const_typeX ;
    typeX x = KokkosArray::create< typeX >( "X" );
    const_typeX xc = x ;

    ASSERT_TRUE( xc == x );
    ASSERT_TRUE( x == xc );

    // typeX xf = xc ; // setting non-const from const must not compile

    check_auto_conversion_to_const( x , x );
  }

  static void run_test_vector()
  {
    enum { Length = 1000 , Count = 8 };

    typedef KokkosArray::View< T[] ,    KokkosArray::LayoutLeft , host > vector_type ;
    typedef KokkosArray::View< T[][0] , KokkosArray::LayoutLeft , host > multivector_type ;

    typedef KokkosArray::View< const T[] ,    KokkosArray::LayoutLeft , host > const_vector_type ;
    typedef KokkosArray::View< const T[][0] , KokkosArray::LayoutLeft , host > const_multivector_type ;

    multivector_type mv = KokkosArray::create<multivector_type>( "mv" , Length , Count );
    vector_type v1 = KokkosArray::view( mv , 0 );
    vector_type v2 = KokkosArray::view( mv , 1 );
    vector_type v3 = KokkosArray::view( mv , 2 );

    ASSERT_TRUE( & v1[0] == & mv(0,0) );
    ASSERT_TRUE( & v2[0] == & mv(0,1) );
    ASSERT_TRUE( & v3[0] == & mv(0,2) );

    const_vector_type cv1( v1 );
    typename vector_type::const_type cv2( v2 );
    typename const_vector_type::const_type ccv2( v2 );

    const_multivector_type cmv( mv );
    typename multivector_type::const_type cmvX( cmv );
    typename const_multivector_type::const_type ccmvX( cmv );
  }
};

} // namespace

/*--------------------------------------------------------------------------*/

