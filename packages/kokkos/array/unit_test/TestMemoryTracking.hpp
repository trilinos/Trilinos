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

class TestMemoryTrackingEntry : public KokkosArray::Impl::MemoryTrackingEntry {
public:

  ~TestMemoryTrackingEntry();

  template< class T >
  TestMemoryTrackingEntry( const T * ptr , const unsigned length , const std::string & label )
    : KokkosArray::Impl::MemoryTrackingEntry( label , typeid(T) , ptr , sizeof(T) * length )
    {}
};

TestMemoryTrackingEntry::~TestMemoryTrackingEntry()
{
  if ( KokkosArray::Impl::MemoryTrackingEntry::count() ) {
    std::cerr << "~TestMemoryTrackingEntry(" ;
    KokkosArray::Impl::MemoryTrackingEntry::print( std::cerr );
    std::cerr << ")" << std::endl ;
  }
}

class TestMemoryTracking {
public:

  TestMemoryTracking() { run_test(); }

  void run_test()
  {
    KokkosArray::Impl::MemoryTracking tracking( "test" );

    int a ;
    double b ;
    long c[10] ;

    tracking.insert( new TestMemoryTrackingEntry( & a , 1 , "a" ) );
    tracking.insert( new TestMemoryTrackingEntry( & b , 1 , "b" ) );
    tracking.insert( new TestMemoryTrackingEntry( c , 10 , "c[10]" ) );

    KokkosArray::Impl::MemoryTrackingEntry * info_a = tracking.query( & a );
    KokkosArray::Impl::MemoryTrackingEntry * info_b = tracking.query( & b );
    KokkosArray::Impl::MemoryTrackingEntry * info_c = tracking.query( c );

    ASSERT_TRUE( 0 != info_a );
    ASSERT_TRUE( 0 != info_b );
    ASSERT_TRUE( 0 != info_c );

    ASSERT_TRUE( info_a->label   == std::string("a") );
    ASSERT_TRUE( info_a->type    == typeid( int ) );
    ASSERT_TRUE( info_a->count() == 1 );
    
    ASSERT_TRUE( info_b->label   == std::string("b") );
    ASSERT_TRUE( info_b->type    == typeid( double ) );
    ASSERT_TRUE( info_b->count() == 1 );
    
    ASSERT_TRUE( info_c->label   == std::string("c[10]") );
    ASSERT_TRUE( info_c->type    == typeid( long ) );
    ASSERT_TRUE( info_c->count() == 1 );

    tracking.increment( & a ); ASSERT_TRUE( 2 == info_a->count() );
    tracking.increment( & a ); ASSERT_TRUE( 3 == info_a->count() );
    tracking.decrement( & a ); ASSERT_TRUE( 2 == info_a->count() );
    tracking.decrement( & a ); ASSERT_TRUE( 1 == info_a->count() );
    tracking.decrement( & a ); ASSERT_TRUE( 0 == tracking.query( & a ) );
    info_a = 0 ;

    tracking.increment( & c[2] ); ASSERT_TRUE( 2 == info_c->count() );
    tracking.decrement( & c[3] ); ASSERT_TRUE( 1 == info_c->count() );
    tracking.decrement( & c[4] ); ASSERT_TRUE( 0 == tracking.query( c ) );
    info_c = 0 ;

    tracking.decrement( & b ); ASSERT_TRUE( 0 == tracking.query( & b ) );
    info_b = 0 ;
  }
};

}

/*--------------------------------------------------------------------------*/



