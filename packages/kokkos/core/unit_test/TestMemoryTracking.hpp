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

#include <gtest/gtest.h>

#include <iostream>
#include <impl/Kokkos_MemoryTracking.hpp>

/*--------------------------------------------------------------------------*/

namespace {

class TestMemoryTracking {
public:

  TestMemoryTracking() { run_test(); }

  void run_test()
  {
    Kokkos::Impl::MemoryTracking<> tracking( "test" );

    int a ;
    double b ;
    long c[10] ;

    tracking.insert( std::string("a") , & a , sizeof(int) , 1 );
    tracking.insert( std::string("b") , & b , sizeof(double) , 1 );
    tracking.insert( std::string("c[10]") , c , sizeof(long) , 10 );

    Kokkos::Impl::MemoryTracking<>::Entry * info_a = tracking.query( & a );
    Kokkos::Impl::MemoryTracking<>::Entry * info_b = tracking.query( & b );
    Kokkos::Impl::MemoryTracking<>::Entry * info_c = tracking.query( c );

    ASSERT_TRUE( 0 != info_a );
    ASSERT_TRUE( 0 != info_b );
    ASSERT_TRUE( 0 != info_c );

    ASSERT_TRUE( info_a->m_alloc_ptr == & a );
    ASSERT_TRUE( info_b->m_alloc_ptr == & b );
    ASSERT_TRUE( info_c->m_alloc_ptr == & c[0] );

    ASSERT_TRUE( info_a->m_type_size == sizeof(int) );
    ASSERT_TRUE( info_b->m_type_size == sizeof(double) );
    ASSERT_TRUE( info_c->m_type_size == sizeof(long) );

    ASSERT_TRUE( info_a->m_array_len == 1 );
    ASSERT_TRUE( info_b->m_array_len == 1 );
    ASSERT_TRUE( info_c->m_array_len == 10 );

    ASSERT_TRUE( std::string( info_a->label() ) == std::string("a") );
    ASSERT_TRUE( std::string( info_b->label() ) == std::string("b") );
    ASSERT_TRUE( std::string( info_c->label() ) == std::string("c[10]") );

    ASSERT_TRUE( info_a->count() == 1 );
    ASSERT_TRUE( info_b->count() == 1 );
    ASSERT_TRUE( info_c->count() == 1 );

    tracking.increment( & a ); ASSERT_TRUE( 2 == info_a->count() );
    tracking.increment( & a ); ASSERT_TRUE( 3 == info_a->count() );

    void * ptr = 0 ;

    ptr = tracking.decrement( & a ); ASSERT_TRUE( ptr == 0 ); ASSERT_TRUE( 2 == info_a->count() );
    ptr = tracking.decrement( & a ); ASSERT_TRUE( ptr == 0 ); ASSERT_TRUE( 1 == info_a->count() );
    ptr = tracking.decrement( & a ); ASSERT_TRUE( ptr == & a ); ASSERT_TRUE( 0 == tracking.query( & a ) );
    info_a = 0 ;

    tracking.increment( & c[2] ); ASSERT_TRUE( 2 == info_c->count() );
    ptr = tracking.decrement( & c[3] ); ASSERT_TRUE( ptr == 0 ); ASSERT_TRUE( 1 == info_c->count() );
    ptr = tracking.decrement( & c[4] ); ASSERT_TRUE( ptr == & c[0] ); ASSERT_TRUE( 0 == tracking.query( c ) );
    info_c = 0 ;

    ptr = tracking.decrement( & b ); ASSERT_TRUE( ptr == & b ); ASSERT_TRUE( 0 == tracking.query( & b ) );
    info_b = 0 ;
  }
};

}

/*--------------------------------------------------------------------------*/



