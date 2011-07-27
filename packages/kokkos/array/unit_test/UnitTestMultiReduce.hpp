/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_ParallelReduce.hpp>

#include <impl/Kokkos_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace UnitTest {

template< typename ScalarType , class DeviceType >
class ReduceMultiFunctorTraits
{
public:

  typedef DeviceType device_type ;

  struct value_type {
    ScalarType value[3] ;
  };

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void init( value_type & dst )
  {
    dst.value[0] = 0 ;
    dst.value[1] = 0 ;
    dst.value[2] = 0 ;
  }

  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  static void join( volatile value_type & dst ,
                    const volatile value_type & src )
  {
    dst.value[0] += src.value[0] ;
    dst.value[1] += src.value[1] ;
    dst.value[2] += src.value[2] ;
  }

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

};


template< typename ScalarType , class DeviceType >
class ReduceMultiFunctor ;

template< typename ScalarType >
class ReduceMultiFunctor< ScalarType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
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

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( size_type iwork , value_type & dst ) const
  {
    const size_type ival = iwork + work_begin ;
    dst.value[0] += 1 ;
    dst.value[1] += 1 + ival ;
    dst.value[2] += work_total - ival ;
  } 

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */

};

} // namespace UnitTest

namespace {

template< typename ScalarType , class DeviceType >
class UnitTestReduceMulti ;

template< typename ScalarType >
class UnitTestReduceMulti< ScalarType , KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE    device_type ;
  typedef device_type::size_type size_type ;

  typedef UnitTest::ReduceMultiFunctorTraits< ScalarType , device_type > reduce_traits ;
  typedef typename reduce_traits::value_type value_type ;

  typedef UnitTest::ReduceMultiFunctor< ScalarType , device_type > functor_type ;

  //------------------------------------

  static std::string name()
  {
    std::string tmp ;
    tmp.append( "UnitTestReduceMulti< Kokkos::" );
    tmp.append( KOKKOS_MACRO_TO_STRING( KOKKOS_MACRO_DEVICE ) );
    tmp.append( " >" );
    return tmp ;
  }

  void error( const char * msg ) const
  {
    std::string tmp = name();
    tmp.append( msg );
    throw std::runtime_error( tmp );
  }

  //------------------------------------

  UnitTestReduceMulti( const size_type nwork ,
                       const size_type nfunctor )
  {
    typedef Kokkos::ValueView< value_type , device_type > result_type ;

    value_type result ;

    Kokkos::MultiFunctorParallelReduce< reduce_traits , result_type >
      reduce_op ;

    for ( size_type j = 0 ; j < nfunctor ; ) {
      const size_type work_beg = (size_t(nwork) * size_t(j) ) / nfunctor ;
      const size_type work_end = (size_t(nwork) * size_t(++j) ) / nfunctor ;
      const size_type work_count = work_end - work_beg ;

      reduce_op.push_back( work_count , functor_type( nwork , work_beg ) );
    }

    reduce_op.result = Kokkos::create_value< result_type >();
    reduce_op.execute();

    Kokkos::deep_copy( result , reduce_op.result );

    const unsigned long nw   = nwork ;
    const unsigned long nsum = nw % 2 ? nw * (( nw + 1 )/2 )
                                      : (nw/2) * ( nw + 1 );

    if ( result.value[0] != (ScalarType) nw ||
         result.value[1] != (ScalarType) nsum ||
         result.value[2] != (ScalarType) nsum ) {
      std::cout << " { " << result.value[0]
                << " , " << result.value[1]
                << " , " << result.value[2]
                << " } != { " << nw
                << " , " << nsum
                << " , " << nsum
                << " }" << std::endl ;
      error( "FAILED" );
    }
  }
};

}

/*--------------------------------------------------------------------------*/

