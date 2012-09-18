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

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <KokkosArray_View.hpp>

/*--------------------------------------------------------------------------*/

namespace {

struct DummyMemorySpace
{
  typedef DummyMemorySpace memory_space ;
  typedef unsigned size_type ;

  static const unsigned align = 64 ; // Byte alignment

  inline static
  size_t preferred_alignment( size_t type_size , size_t n )
  {
    // goal: ( n * type_size ) % align == 0
    while ( 0 != ( n * type_size ) % align ) ++n ;
    return n ;
  }
};

/*--------------------------------------------------------------------------*/

template< class Layout , class Type >
struct DefineShape {
  typedef typename KokkosArray::Impl::AnalyzeShape<Type,Layout>::shape type ;
};

template< class Type >
struct ExtractValueType {
  typedef typename KokkosArray::Impl::AnalyzeShape<Type>::value_type type ;
};

template< class Type >
struct ArrayType { typedef Type type ; };

template < class Device >
void test_view_impl()
{
  typedef typename Device::memory_space memory_space ;

  typedef ArrayType< int[100]                >::type type_01 ;
  typedef ArrayType< int[]                   >::type type_11 ;
  typedef ArrayType< int[5][6][700]          >::type type_03 ;
  typedef ArrayType< double[][8][9][900]     >::type type_14 ;
  typedef ArrayType< long**                  >::type type_22 ;
  typedef ArrayType< short***[5][6][7]       >::type type_36 ;
  typedef ArrayType< const short***[5][6][7] >::type const_type_36 ;
  typedef ArrayType< short **[5][6][7]       >::type type_25 ;
  typedef ArrayType< const short **[5][6][7] >::type const_type_25 ;

  typedef typename KokkosArray::Impl::StaticAssertSame<
     typename KokkosArray::Impl::AnalyzeShape<type_25>::const_type ,
     typename KokkosArray::Impl::AnalyzeShape<const_type_25>::type 
       > ok_const_25 ;

  ASSERT_TRUE( ( KokkosArray::Impl::is_same< ExtractValueType<type_03>::type , int >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::is_same< ExtractValueType<type_14>::type , double >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::is_same< ExtractValueType<type_22>::type , long >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::is_same< ExtractValueType<type_36>::type , short >::value ) );

  ASSERT_FALSE( ( KokkosArray::Impl::is_same< ExtractValueType<type_36>::type , int >::value ) );


  typedef typename DefineShape< KokkosArray::LayoutLeft , type_01>::type  shape_01_type ;
  typedef typename DefineShape< KokkosArray::LayoutLeft , type_11>::type  shape_11_type ;
  typedef typename DefineShape< KokkosArray::LayoutLeft , type_03>::type  shape_03_type ;
  typedef typename DefineShape< KokkosArray::LayoutRight, type_14>::type  shape_14_type ;
  typedef typename DefineShape< KokkosArray::LayoutRight, type_22>::type  shape_22_type ;
  typedef typename DefineShape< KokkosArray::LayoutRight, type_36>::type  shape_36_type ;

  ASSERT_TRUE( ( KokkosArray::Impl::StaticAssert< shape_36_type::rank == 6 >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::StaticAssert< shape_03_type::rank == 3 >::value ) );

  shape_01_type shape_01 = shape_01_type::create<DummyMemorySpace>();
  shape_11_type shape_11 = shape_11_type::create<DummyMemorySpace>( 1000 );
  shape_03_type shape_03 = shape_03_type::create<DummyMemorySpace>();
  shape_14_type shape_14 = shape_14_type::create<DummyMemorySpace>( 0 );
  shape_22_type shape_22 = shape_22_type::create<DummyMemorySpace>( 0 , 0 );
  shape_36_type shape_36 = shape_36_type::create<DummyMemorySpace>( 10 , 20 , 30 );

  ASSERT_TRUE( shape_01.rank_dynamic == 0u );
  ASSERT_TRUE( shape_01.rank         == 1u );
  ASSERT_TRUE( shape_01.N0           == 100u );

  ASSERT_TRUE( shape_11.rank_dynamic == 1u );
  ASSERT_TRUE( shape_11.rank         == 1u );
  ASSERT_TRUE( shape_11.N0           == 1000u );

  ASSERT_TRUE( shape_03.rank_dynamic == 0u );
  ASSERT_TRUE( shape_03.rank         == 3u );
  ASSERT_TRUE( shape_03.N0           == 5u );
  ASSERT_TRUE( shape_03.N1           == 6u );
  ASSERT_TRUE( shape_03.N2           == 700u );

  ASSERT_TRUE( shape_14.rank_dynamic == 1u );
  ASSERT_TRUE( shape_14.rank         == 4u );
  ASSERT_TRUE( shape_14.StaticN0     == 0u );
  ASSERT_TRUE( shape_14.N0           == 0u );
  ASSERT_TRUE( shape_14.N1           == 8u );
  ASSERT_TRUE( shape_14.N2           == 9u );
  ASSERT_TRUE( shape_14.N3           == 900u );

  ASSERT_TRUE( shape_22.rank_dynamic == 2u );
  ASSERT_TRUE( shape_22.rank         == 2u );
  ASSERT_TRUE( shape_22.StaticN0     == 0u );
  ASSERT_TRUE( shape_22.StaticN1     == 0u );
  ASSERT_TRUE( shape_22.N0           == 0u );
  ASSERT_TRUE( shape_22.N1           == 0u );

  ASSERT_TRUE( shape_36.rank_dynamic == 3u );
  ASSERT_TRUE( shape_36.rank         == 6u );
  ASSERT_TRUE( shape_36.StaticN0     == 0u );
  ASSERT_TRUE( shape_36.StaticN1     == 0u );
  ASSERT_TRUE( shape_36.StaticN2     == 0u );
  ASSERT_TRUE( shape_36.N0           == 10u );
  ASSERT_TRUE( shape_36.N1           == 20u );
  ASSERT_TRUE( shape_36.N2           == 30u );
  ASSERT_TRUE( shape_36.N3           == 5u  );
  ASSERT_TRUE( shape_36.N4           == 6u  );
  ASSERT_TRUE( shape_36.N5           == 7u  );


  ASSERT_TRUE( shape_01 == shape_01 );
  ASSERT_TRUE( shape_11 == shape_11 );
  ASSERT_TRUE( shape_36 == shape_36 );
  ASSERT_TRUE( shape_01 != shape_36 );
  ASSERT_TRUE( shape_22 != shape_36 );

  //------------------------------------------------------------------------
}

} /* namespace */

/*--------------------------------------------------------------------------*/

