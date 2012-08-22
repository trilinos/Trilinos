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

#ifndef KOKKOSARRAY_MACRO_DEVICE
#error "KOKKOSARRAY_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <KokkosArray_ParallelReduce.hpp>

/*--------------------------------------------------------------------------*/

namespace Test {

template< typename ScalarType , class DeviceType >
class ReduceMultiFunctorTraits
{
public:

  typedef DeviceType device_type ;

  struct value_type {
    ScalarType value[3] ;
  };

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  static void init( value_type & dst )
  {
    dst.value[0] = 0 ;
    dst.value[1] = 0 ;
    dst.value[2] = 0 ;
  }

  inline
  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & dst ,
                    const volatile value_type & src )
  {
    dst.value[0] += src.value[0] ;
    dst.value[1] += src.value[1] ;
    dst.value[2] += src.value[2] ;
  }

#endif /* defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION ) */

};


template< typename ScalarType , class DeviceType >
class ReduceMultiFunctor ;

template< typename ScalarType >
class ReduceMultiFunctor< ScalarType , KOKKOSARRAY_MACRO_DEVICE >
{
public:
  typedef KOKKOSARRAY_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  typedef ReduceMultiFunctorTraits< ScalarType , device_type > reduce_traits ;
  typedef typename reduce_traits::value_type value_type ;

  const size_type work_total ;
  const size_type work_begin ;

  ReduceMultiFunctor( const size_type & arg_total ,
                      const size_type & arg_begin )
    : work_total( arg_total )
    , work_begin( arg_begin )
    {}

  ReduceMultiFunctor( const ReduceMultiFunctor & rhs )
    : work_total( rhs.work_total )
    , work_begin( rhs.work_begin )
    {}

#if defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION )

  KOKKOSARRAY_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork , value_type & dst ) const
  {
    const size_type ival = iwork + work_begin ;
    dst.value[0] += 1 ;
    dst.value[1] += 1 + ival ;
    dst.value[2] += work_total - ival ;
  }

#endif /* defined( KOKKOSARRAY_MACRO_DEVICE_FUNCTION ) */

};

} // namespace Test

namespace {

template< typename ScalarType , class DeviceType >
class TestReduceMulti ;

template< typename ScalarType >
class TestReduceMulti< ScalarType , KOKKOSARRAY_MACRO_DEVICE >
{
public:
  typedef KOKKOSARRAY_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  typedef Test::ReduceMultiFunctorTraits< ScalarType , device_type > reduce_traits ;
  typedef typename reduce_traits::value_type value_type ;

  typedef Test::ReduceMultiFunctor< ScalarType , device_type > functor_type ;


  //------------------------------------
  TestReduceMulti( const size_type nwork ,
                       const size_type nfunctor )
  { run_test(nwork, nfunctor); }

  void run_test( const size_type nwork ,
                       const size_type nfunctor )
  {
    value_type result ;

    { // Destruction of the 'result_functor_type' copies result
      // data from the device to the host, as necessary.

      typedef KokkosArray::Impl
                ::FunctorAssignment< value_type , device_type >
                  result_functor_type ;

      result_functor_type result_functor( result );

      KokkosArray::MultiFunctorParallelReduce< reduce_traits , result_functor_type , device_type >
        reduce_op( result_functor );

      for ( size_type j = 0 ; j < nfunctor ; ) {
        const size_type work_beg = (size_t(nwork) * size_t(j) ) / nfunctor ;
        const size_type work_end = (size_t(nwork) * size_t(++j) ) / nfunctor ;
        const size_type work_count = work_end - work_beg ;

        reduce_op.push_back( work_count , functor_type( nwork , work_beg ) );
      }

      reduce_op.execute();
    }

    const unsigned long nw   = nwork ;
    const unsigned long nsum = nw % 2 ? nw * (( nw + 1 )/2 )
                                      : (nw/2) * ( nw + 1 );

    ASSERT_EQ( result.value[0], (ScalarType) nw);
    ASSERT_EQ( result.value[1], (ScalarType) nsum);
    ASSERT_EQ( result.value[2], (ScalarType) nsum);
  }
};

}

/*--------------------------------------------------------------------------*/

