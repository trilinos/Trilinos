/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <KokkosArray_ParallelReduce.hpp>

/*--------------------------------------------------------------------------*/

namespace Test {

template< typename ScalarType , class DeviceType >
class ReduceFunctor
{
public:
  typedef DeviceType  device_type ;
  typedef typename device_type::size_type size_type ;

  struct value_type {
    ScalarType value[3] ;
  };

  const size_type nwork ;

  ReduceFunctor( const size_type & arg_nwork ) : nwork( arg_nwork ) {}

  ReduceFunctor( const ReduceFunctor & rhs )
    : nwork( rhs.nwork ) {}

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  static void init( value_type & dst )
  {
    dst.value[0] = 0 ;
    dst.value[1] = 0 ;
    dst.value[2] = 0 ;
  }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  static void join( volatile value_type & dst ,
                    const volatile value_type & src )
  {
    dst.value[0] += src.value[0] ;
    dst.value[1] += src.value[1] ;
    dst.value[2] += src.value[2] ;
  }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void operator()( size_type iwork , value_type & dst ) const
  {
    dst.value[0] += 1 ;
    dst.value[1] += iwork + 1 ;
    dst.value[2] += nwork - iwork ;
  }

};

template< typename ScalarType , class DeviceType >
class RuntimeReduceFunctor
{
public:
  // Required for functor:
  typedef DeviceType  device_type ;
  typedef ScalarType  value_type[] ;
  const unsigned      value_count ;


  // Unit test details:

  typedef typename device_type::size_type  size_type ;

  const size_type     nwork ;

  RuntimeReduceFunctor( const size_type & arg_nwork )
    : value_count( 3 )
    , nwork( arg_nwork ) {}

  RuntimeReduceFunctor( const RuntimeReduceFunctor & rhs )
    : value_count( rhs.value_count )
    , nwork( rhs.nwork ) {}

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  static void init( value_type dst ,
                    const size_type )
  {
    dst[0] = 0 ;
    dst[1] = 0 ;
    dst[2] = 0 ;
  }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  static void join( volatile ScalarType dst[] ,
                    const volatile ScalarType src[] ,
                    const size_type n )
  {
    dst[0] += src[0] ;
    dst[1] += src[1] ;
    dst[2] += src[2] ;
  }

  KOKKOSARRAY_INLINE_DEVICE_FUNCTION
  void operator()( size_type iwork , ScalarType dst[] ) const
  {
    dst[0] += 1 ;
    dst[1] += iwork + 1 ;
    dst[2] += nwork - iwork ;
  }
};

} // namespace Test

namespace {

template< typename ScalarType , class DeviceType >
class TestReduce
{
public:
  typedef DeviceType    device_type ;
  typedef typename device_type::size_type size_type ;

  //------------------------------------

  TestReduce( const size_type & nwork )
  {
    run_test(nwork);
  }

  void run_test( const size_type & nwork )
  {
    typedef Test::ReduceFunctor< ScalarType , device_type > functor_type ;
    typedef typename functor_type::value_type value_type ;

    enum { Repeat = 100 };

    value_type result[ Repeat ];

    const unsigned long nw   = nwork ;
    const unsigned long nsum = nw % 2 ? nw * (( nw + 1 )/2 )
                                      : (nw/2) * ( nw + 1 );

    for ( unsigned i = 0 ; i < Repeat ; ++i ) {
      result[i] = KokkosArray::parallel_reduce( nwork , functor_type(nwork) );
    }

    for ( unsigned i = 0 ; i < Repeat ; ++i ) {
      ASSERT_EQ( result[i].value[0], (ScalarType) nw);
      ASSERT_EQ( result[i].value[1], (ScalarType) nsum);
      ASSERT_EQ( result[i].value[2], (ScalarType) nsum);
    }
  }
};

template< typename ScalarType , class DeviceType >
class TestReduceDynamic
{
public:
  typedef DeviceType    device_type ;
  typedef typename device_type::size_type size_type ;

  //------------------------------------

  TestReduceDynamic( const size_type & nwork )
  {
    run_test_dynamic(nwork);
  }

  void run_test_dynamic( const size_type & nwork )
  {
    typedef Test::RuntimeReduceFunctor< ScalarType , device_type > functor_type ;
    typedef KokkosArray::Impl::ParallelReduceFunctorValue< ScalarType[] , device_type > finalize_type ;

    enum { Repeat = 100 };

    ScalarType result[ Repeat ][3] ;

    const finalize_type finalize(3);

    const unsigned long nw   = nwork ;
    const unsigned long nsum = nw % 2 ? nw * (( nw + 1 )/2 )
                                      : (nw/2) * ( nw + 1 );

    for ( unsigned i = 0 ; i < Repeat ; ++i ) {

//      KokkosArray::parallel_reduce( nwork , functor_type(nwork) , result[i]  );

      KokkosArray::parallel_reduce( nwork , functor_type(nwork) , finalize );

      finalize.result( result[i] );
    }

    for ( unsigned i = 0 ; i < Repeat ; ++i ) {
      ASSERT_EQ( result[i][0], (ScalarType) nw);
      ASSERT_EQ( result[i][1], (ScalarType) nsum);
      ASSERT_EQ( result[i][2], (ScalarType) nsum);
    }
  }
};

}

/*--------------------------------------------------------------------------*/

