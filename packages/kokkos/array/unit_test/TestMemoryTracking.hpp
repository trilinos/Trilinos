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

#include <gtest/gtest.h>

#include <iostream>
#include <impl/KokkosArray_MemoryTracking.hpp>

/*--------------------------------------------------------------------------*/

namespace {

class TestMemoryTracking {
public:

  TestMemoryTracking() { run_test(); }

  void run_test()
  {
    KokkosArray::Impl::MemoryTracking tracking( "test" );
    KokkosArray::Impl::MemoryTracking::Info info ;

    int a ;
    double b ;
    long c[10] ;

    tracking.track( & a , 1 , "a" );
    tracking.track( & b , 1 , "b" );
    tracking.track( c , 10 , "c[10]" );

    info = tracking.query( & a );
    ASSERT_TRUE( info.label  == std::string("a") );
    ASSERT_TRUE( info.type   == & typeid( int ) );
    ASSERT_TRUE( info.size   == sizeof( int ) );
    ASSERT_TRUE( info.length == 1 );
    ASSERT_TRUE( info.count  == 1 );
    
    info = tracking.query( & b );
    ASSERT_TRUE( info.label  == std::string("b") );
    ASSERT_TRUE( info.type   == & typeid( double ) );
    ASSERT_TRUE( info.size   == sizeof( double ) );
    ASSERT_TRUE( info.length == 1 );
    ASSERT_TRUE( info.count  == 1 );
    
    info = tracking.query( c );
    ASSERT_TRUE( info.label  == std::string("c[10]") );
    ASSERT_TRUE( info.type   == & typeid( long ) );
    ASSERT_TRUE( info.size   == sizeof( long ) );
    ASSERT_TRUE( info.length == 10 );
    ASSERT_TRUE( info.count  == 1 );

    tracking.increment( & a ); ASSERT_TRUE( 2 == tracking.query( & a ).count );
    tracking.increment( & a ); ASSERT_TRUE( 3 == tracking.query( & a ).count );

    info = tracking.query( & a );
    ASSERT_TRUE( info.label  == std::string("a") );
    ASSERT_TRUE( info.type   == & typeid( int ) );
    ASSERT_TRUE( info.size   == sizeof( int ) );
    ASSERT_TRUE( info.length == 1 );
    ASSERT_TRUE( info.count  == 3 );
    
    info = tracking.query( & b );
    ASSERT_TRUE( info.label  == std::string("b") );
    ASSERT_TRUE( info.type   == & typeid( double ) );
    ASSERT_TRUE( info.size   == sizeof( double ) );
    ASSERT_TRUE( info.length == 1 );
    ASSERT_TRUE( info.count  == 1 );
    
    ASSERT_TRUE( 0 != tracking.decrement( & b ) );

    info = tracking.query( & b );
    ASSERT_TRUE( info.label  == std::string() );
    ASSERT_TRUE( info.begin  == 0 );
    ASSERT_TRUE( info.type   == 0 );
    ASSERT_TRUE( info.size   == 0 );
    ASSERT_TRUE( info.length == 0 );
    ASSERT_TRUE( info.count  == 0 );

    ASSERT_TRUE( 0 == tracking.decrement( & a ) );
    ASSERT_TRUE( 2 == tracking.query( & a ).count );

    info = tracking.query( & a );
    ASSERT_TRUE( info.label  == std::string("a") );
    ASSERT_TRUE( info.type   == & typeid( int ) );
    ASSERT_TRUE( info.size   == sizeof( int ) );
    ASSERT_TRUE( info.length == 1 );
    ASSERT_TRUE( info.count  == 2 );
    
    ASSERT_TRUE( 0 == tracking.decrement( & a ) );
    ASSERT_TRUE( 1 == tracking.query( & a ).count );
    ASSERT_TRUE( 0 != tracking.decrement( & a ) );

    info = tracking.query( & a );
    ASSERT_TRUE( info.label  == std::string() );
    ASSERT_TRUE( info.begin  == 0 );
    ASSERT_TRUE( info.type   == 0 );
    ASSERT_TRUE( info.size   == 0 );
    ASSERT_TRUE( info.length == 0 );
    ASSERT_TRUE( info.count  == 0 );

    tracking.increment( & c[2] ); ASSERT_EQ( 2u , tracking.query( & c[2] ).count );
    ASSERT_TRUE( 0 == tracking.decrement( & c[3] ) );
    ASSERT_EQ( 1u , tracking.query( & c[3] ).count );

    ASSERT_TRUE( 0 != tracking.decrement( c ) );
    ASSERT_TRUE( tracking.empty() );
  }
};

}

/*--------------------------------------------------------------------------*/



