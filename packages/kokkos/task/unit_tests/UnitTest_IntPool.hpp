/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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


#ifndef KOKKOS_UNIT_TEST_INTPOOL_HPP
#define KOKKOS_UNIT_TEST_INTPOOL_HPP

#include <iostream>
#include <impl/Kokkos_IntPool.hpp>

namespace Test {

//----------------------------------------------------------------------------

template< class Device >
struct TestIntPool {

  typedef Device device_type ;

  typedef Kokkos::Impl::IntPool<device_type> IntPool ;


  KOKKOS_INLINE_FUNCTION
  void operator()( int i ) const
    {
      typename IntPool::ClaimResult result ;

      result.value = i ;

      do {
        result = x.claim( result.value );
      } while( result.status == IntPool::ClaimResult::FAIL );
    }

  IntPool x ;
};


//----------------------------------------------------------------------------

template< class Device >
void test_intpool( const int );

template<>
void test_intpool<Kokkos::Serial>( const int n )
{
  typedef Kokkos::Impl::IntPool<Kokkos::Serial> IntPool ;

  {
    IntPool x( n );

    for ( int i = 0 ; i < n ; ++i ) {

      int value = i ;

      IntPool::ClaimStatus result = x.claim(value);

      if ( IntPool::SUCCESS != result || i != value ) {
        std::cerr << "test claim FAILED " << i << std::endl ;
        x.print( std::cerr );
        return ;
      }

      x.release( value );

      result = x.claim( value );

      if ( IntPool::SUCCESS != result || i != value ) {
        std::cerr << "test claim FAILED " << i << std::endl ;
        x.print( std::cerr );
        return ;
      }

      x[ value ] = 10*(i+1);
    }

    int value = 0 ;
    if ( x.claim( value ) != IntPool::FULL ) {
      std::cerr << "test claim full FAILED " << std::endl ;
      x.print( std::cerr );
      return ;
    }
  }

  {
    IntPool x( n );

    int value = 0 ;
    for ( int i = 0 ; i < n ; ++i ) {
      x.claim( value );
      x[ value ] = 10 * ( value + 1 );
    }
    for ( int i = 0 ; i < n ; i += 2 ) { x.release( i ); }

    value = 0 ;
    for ( int i = 0 ; i < n ; i += 2 ) {
      if ( IntPool::SUCCESS != x.claim( value ) ) {
        std::cerr << "test claim repeat failed " << i << std::endl ;
        x.print( std::cerr );
        return ;
      }
      x[ value ] = 10 * ( value + 1 );
    }
    for ( int i = 0 ; i < n ; ++i ) {
      if ( ! x.is_claimed(i) || x[i] != 10*(i+1) ) {
        std::cerr << "test claim repeat value failed " << i << std::endl ;
        x.print( std::cerr );
        return ;
      }
    }
  }

  {
    IntPool x( n );

    int value = 0 ;
    for ( int i = 0 ; i < n ; i += 2 ) {

      IntPool::ClaimStatus result = x.claim( value );
    
      if ( IntPool::SUCCESS != result || i != value ) {
        std::cerr << "test claim FAILED " << i << std::endl ;
        x.print( std::cerr );
        return ;
      }

      if ( i + 1 < n ) {
        value = i ;
        result = x.claim( value );

        if ( IntPool::SUCCESS != result || i + 1 != value ) {
          std::cerr << "test claim FAILED " << i + 1 << " " << value << std::endl ;
          x.print( std::cerr );
          return ;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------

} // namespace Test

#endif /* #ifndef KOKKOS_UNIT_TEST_INTPOOL_HPP */

