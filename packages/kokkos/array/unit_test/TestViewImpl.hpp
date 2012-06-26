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

template < class ShapeType ,
           class MemorySpace ,
           class Layout = typename ShapeType::array_layout ,
           unsigned Rank = ShapeType::rank >
struct TestShapeMap ;

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutRight , 8 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  for ( unsigned i5 = 0 ; i5 < shape.N5 ; ++i5 )
  for ( unsigned i6 = 0 ; i6 < shape.N6 ; ++i6 )
  for ( unsigned i7 = 0 ; i7 < shape.N7 ; ++i7 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4, i5, i6, i7 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutRight , 7 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  for ( unsigned i5 = 0 ; i5 < shape.N5 ; ++i5 )
  for ( unsigned i6 = 0 ; i6 < shape.N6 ; ++i6 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4, i5, i6 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutRight , 6 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  for ( unsigned i5 = 0 ; i5 < shape.N5 ; ++i5 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4, i5 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};


template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutRight , 5 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutRight , 4 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutRight , 3 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutRight , 2 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  {
    const size_t j = shape_map::offset( shape, i0, i1 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

/*--------------------------------------------------------------------------*/

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutLeft , 8 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i7 = 0 ; i7 < shape.N7 ; ++i7 )
  for ( unsigned i6 = 0 ; i6 < shape.N6 ; ++i6 )
  for ( unsigned i5 = 0 ; i5 < shape.N5 ; ++i5 )
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4, i5, i6, i7 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutLeft , 7 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i6 = 0 ; i6 < shape.N6 ; ++i6 )
  for ( unsigned i5 = 0 ; i5 < shape.N5 ; ++i5 )
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4, i5, i6 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutLeft , 6 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i5 = 0 ; i5 < shape.N5 ; ++i5 )
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4, i5 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutLeft , 5 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i4 = 0 ; i4 < shape.N4 ; ++i4 )
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3, i4 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutLeft , 4 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i3 = 0 ; i3 < shape.N3 ; ++i3 )
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2, i3 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutLeft , 3 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i2 = 0 ; i2 < shape.N2 ; ++i2 )
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  {
    const size_t j = shape_map::offset( shape, i0, i1, i2 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

template < class ShapeType ,
           class MemorySpace >
struct TestShapeMap< ShapeType , MemorySpace , KokkosArray::LayoutLeft , 2 > {

static void apply( const ShapeType & shape )
{
  typedef KokkosArray::Impl::ShapeMap< ShapeType , MemorySpace > shape_map ;

  const size_t alloc = KokkosArray::Impl::allocation_count( shape );

  ASSERT_TRUE( KokkosArray::Impl::cardinality_count( shape ) <= alloc );

  long offset = -1 ;
  for ( unsigned i1 = 0 ; i1 < shape.N1 ; ++i1 ) 
  for ( unsigned i0 = 0 ; i0 < shape.N0 ; ++i0 )
  {
    const size_t j = shape_map::offset( shape, i0, i1 );
    ASSERT_TRUE( j < alloc );
    ASSERT_TRUE( offset < (long) j );
    offset = j ;
  }
}
};

/*--------------------------------------------------------------------------*/

template< class Type , class MemorySpace >
struct TestLeftAndRight
{
  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutRight ,
      typename KokkosArray::Impl::change_empty_extent_to_zero_extent<Type>::type
    > right_type ;
  
  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutLeft ,
      typename KokkosArray::Impl::change_empty_extent_to_zero_extent<Type>::type
    > left_type ;

  static void apply()
  {
    left_type left = KokkosArray::Impl::Factory<left_type,DummyMemorySpace>::create();
    right_type right = KokkosArray::Impl::Factory<right_type,DummyMemorySpace>::create();

    TestShapeMap< left_type ,  MemorySpace >::apply( left );
    TestShapeMap< right_type , MemorySpace >::apply( right );
  }
};

/*--------------------------------------------------------------------------*/

template < class Device >
void test_view_impl()
{
  typedef typename Device::memory_space_new memory_space ;

  typedef int    type_01 [100];
  typedef int    type_11 [0];
  typedef int    type_03 [5][6][700];
  typedef double type_14 [0][8][9][900];
  typedef long   type_22 [0][0] ;
  typedef short  type_36 [][0][0][5][6][7] ;
  typedef short  type_25 [0][0][5][6][7] ;

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<int[]>::value == 1u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        int[]>::value == 1u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<int[][0][5]>::value == 2u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        int[][0][5]>::value == 3u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<int[0]>::value == 1u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        int[0]>::value == 1u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<int[0][0][5]>::value == 2u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        int[0][0][5]>::value == 3u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<type_01>::value == 0u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        type_01>::value == 1u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<type_11>::value == 1u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        type_11>::value == 1u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<type_03>::value == 0u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        type_03>::value == 3u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<type_14>::value == 1u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        type_14>::value == 4u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<type_22>::value == 2u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        type_22>::value == 2u >::value );

  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank_dynamic<type_36>::value == 3u >::value );
  ASSERT_TRUE( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::rank<        type_36>::value == 6u >::value );

  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,0>::value == 5u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,1>::value == 6u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,2>::value == 700u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,3>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,4>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,5>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,6>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,7>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_03,8>::value == 0u ) );

  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,0>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,1>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,2>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,3>::value == 5u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,4>::value == 6u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,5>::value == 7u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,6>::value == 0u ) );
  ASSERT_TRUE( ( KokkosArray::Impl::extent<type_36,7>::value == 0u ) );

  ASSERT_TRUE( ( KokkosArray::Impl::is_same< KokkosArray::Impl::remove_all_extents<type_03>::type , int >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::is_same< KokkosArray::Impl::remove_all_extents<type_14>::type , double >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::is_same< KokkosArray::Impl::remove_all_extents<type_22>::type , long >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::is_same< KokkosArray::Impl::remove_all_extents<type_36>::type , short >::value ) );

  ASSERT_FALSE( ( KokkosArray::Impl::is_same< KokkosArray::Impl::remove_all_extents<type_36>::type , int >::value ) );

  ASSERT_TRUE( ( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::Shape< KokkosArray::LayoutLeft , type_36 >::rank == 6 >::value ) );
  ASSERT_TRUE( ( KokkosArray::Impl::StaticAssert< KokkosArray::Impl::Shape< KokkosArray::LayoutRight, type_36 >::rank == 6 >::value ) );

  typedef typename KokkosArray::Impl::StaticAssertSame< typename KokkosArray::Impl::remove_extent< type_36 , 0 >::type ,
                                                        type_25 >::type ok_remove_extent_0 ;

  typedef typename KokkosArray::Impl::StaticAssertSame< typename KokkosArray::Impl::remove_extent< type_36 , 2 >::type ,
                                                        short[0][0] [5][6][7] >::type ok_remove_extent_2 ;

  typedef typename KokkosArray::Impl::StaticAssertSame< typename KokkosArray::Impl::remove_extent< type_36 , 3 >::type ,
                                                        short[0][0][0] [6][7] >::type ok_remove_extent_3 ;

  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutLeft ,
      KokkosArray::Impl::change_empty_extent_to_zero_extent<type_01>::type
    > shape_01_type ;

  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutLeft ,
      KokkosArray::Impl::change_empty_extent_to_zero_extent<type_11>::type
    > shape_11_type ;

  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutLeft ,
      KokkosArray::Impl::change_empty_extent_to_zero_extent<type_03>::type
    > shape_03_type ;

  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutRight ,
      KokkosArray::Impl::change_empty_extent_to_zero_extent<type_14>::type
    > shape_14_type ;

  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutRight ,
      KokkosArray::Impl::change_empty_extent_to_zero_extent<type_22>::type
    > shape_22_type ;

  typedef KokkosArray::Impl::Shape<
      KokkosArray::LayoutRight ,
      KokkosArray::Impl::change_empty_extent_to_zero_extent<type_36>::type
    > shape_36_type ;

  shape_01_type shape_01 = KokkosArray::Impl::Factory<shape_01_type,DummyMemorySpace>::create();
  shape_11_type shape_11 = KokkosArray::Impl::Factory<shape_11_type,DummyMemorySpace>::create( 1000 );
  shape_03_type shape_03 = KokkosArray::Impl::Factory<shape_03_type,DummyMemorySpace>::create();
  shape_14_type shape_14 = KokkosArray::Impl::Factory<shape_14_type,DummyMemorySpace>::create( 0 );
  shape_22_type shape_22 = KokkosArray::Impl::Factory<shape_22_type,DummyMemorySpace>::create( 0 , 0 );
  shape_36_type shape_36 = KokkosArray::Impl::Factory<shape_36_type,DummyMemorySpace>::create( 10 , 20 , 30 );

  ASSERT_TRUE( shape_01.rank_dynamic == 0u );
  ASSERT_TRUE( shape_01.rank         == 1u );

  ASSERT_TRUE( shape_11.rank_dynamic == 1u );
  ASSERT_TRUE( shape_11.rank         == 1u );

  ASSERT_TRUE( shape_03.rank_dynamic == 0u );
  ASSERT_TRUE( shape_03.rank         == 3u );

  ASSERT_TRUE( shape_14.rank_dynamic == 1u );
  ASSERT_TRUE( shape_14.rank         == 4u );

  ASSERT_TRUE( shape_22.rank_dynamic == 2u );
  ASSERT_TRUE( shape_22.rank         == 2u );

  ASSERT_TRUE( shape_36.rank_dynamic == 3u );
  ASSERT_TRUE( shape_36.rank         == 6u );

  ASSERT_TRUE( shape_01.N0 == 100u );
  ASSERT_TRUE( shape_11.N0 == 1000u );

  ASSERT_TRUE( shape_03.N0 == 5u );
  ASSERT_TRUE( shape_03.N1 == 6u );
  ASSERT_TRUE( shape_03.N2 == 700u );
  ASSERT_TRUE( shape_03.N3 == 0u );
  ASSERT_TRUE( shape_03.N4 == 0u );
  ASSERT_TRUE( shape_03.N5 == 0u );
  ASSERT_TRUE( shape_03.N6 == 0u );
  ASSERT_TRUE( shape_03.N7 == 0u );

  ASSERT_TRUE( shape_36.N0 == 10u );
  ASSERT_TRUE( shape_36.N1 == 20u );
  ASSERT_TRUE( shape_36.N2 == 30u );
  ASSERT_TRUE( shape_36.N3 == 5u  );
  ASSERT_TRUE( shape_36.N4 == 6u  );
  ASSERT_TRUE( shape_36.N5 == 7u  );
  ASSERT_TRUE( shape_36.N6 == 0u  );
  ASSERT_TRUE( shape_36.N7 == 0u  );

  ASSERT_TRUE( shape_01 == shape_01 );
  ASSERT_TRUE( shape_11 == shape_11 );
  ASSERT_TRUE( shape_36 == shape_36 );
  ASSERT_TRUE( shape_01 != shape_36 );
  ASSERT_TRUE( shape_22 != shape_36 );

  TestShapeMap< shape_03_type , memory_space >::apply( shape_03 );
  TestShapeMap< shape_36_type , memory_space >::apply( shape_36 );

  TestLeftAndRight< int[2][3] , memory_space >::apply();
  TestLeftAndRight< int[2][3][4] , memory_space >::apply();
  TestLeftAndRight< int[2][3][4][2] , memory_space >::apply();
  TestLeftAndRight< int[2][3][4][2][3] , memory_space >::apply();
  TestLeftAndRight< int[2][3][4][2][3][4] , memory_space >::apply();
  TestLeftAndRight< int[2][3][4][2][3][4][2] , memory_space >::apply();
  TestLeftAndRight< int[2][3][4][2][3][4][2][3] , memory_space >::apply();

  // KokkosArray::View< int , KokkosArray::LayoutLeft , DummyMemorySpace > compile_error ;
}

} /* namespace */

/*--------------------------------------------------------------------------*/

