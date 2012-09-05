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

#ifndef KOKKOSARRAY_SHAPE_HPP
#define KOKKOSARRAY_SHAPE_HPP

#include <typeinfo>
#include <KokkosArray_Layout.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  The shape of a Kokkos Array with dynamic and static dimensions.
 *          Dynamic dimensions are member values and static dimensions are
 *          'static const' values.
 *
 *  The upper bound on the array rank is eight.
 */
template< class Layout ,
          unsigned ValueSize ,
          unsigned Rank ,
          unsigned s0  = 1 ,
          unsigned s1  = 1 ,
          unsigned s2  = 1 ,
          unsigned s3  = 1 ,
          unsigned s4  = 1 ,
          unsigned s5  = 1 ,
          unsigned s6  = 1 ,
          unsigned s7  = 1 >
struct Shape ;

template< class ShapeType >
struct ShapeMap ;

//----------------------------------------------------------------------------
/** \brief  Shape equality if the value type, layout, and dimensions
 *          are equal.
 */
template< class    xLayout , unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          class    yLayout , unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
bool operator == ( const Shape<xLayout,xSize,xRank,
                               xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
                   const Shape<yLayout,ySize,yRank,
                               yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  enum { same_size = xSize == ySize };
  enum { same_rank = xRank == yRank };

  // the array layout only matters for 1 < rank
  enum { same_layout = xRank < 2 || is_same< xLayout , yLayout >::value };

  return same_size && same_layout && same_rank &&
         x.N0 == y.N0 && x.N1 == y.N1 && x.N2 == y.N2 && x.N3 == y.N3 &&
         x.N4 == y.N4 && x.N5 == y.N5 && x.N6 == y.N6 && x.N7 == y.N7 &&
         x.Stride == y.Stride ;
}

template< class    xLayout , unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          class    yLayout , unsigned ySize ,unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
bool operator != ( const Shape<xLayout,xSize,xRank,
                               xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
                   const Shape<yLayout,ySize,yRank,
                               yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{ return ! operator == ( x , y ); }

//----------------------------------------------------------------------------

void assert_shapes_are_equal_throw(
  const std::type_info & x_layout ,
  const unsigned x_value_size ,
  const unsigned x_rank , const unsigned x_stride ,
  const unsigned x_N0 , const unsigned x_N1 ,
  const unsigned x_N2 , const unsigned x_N3 ,
  const unsigned x_N4 , const unsigned x_N5 ,
  const unsigned x_N6 , const unsigned x_N7 ,

  const std::type_info & y_layout ,
  const unsigned y_value_size ,
  const unsigned y_rank , const unsigned y_stride ,
  const unsigned y_N0 , const unsigned y_N1 ,
  const unsigned y_N2 , const unsigned y_N3 ,
  const unsigned y_N4 , const unsigned y_N5 ,
  const unsigned y_N6 , const unsigned y_N7 );

template< class    xLayout , unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          class    yLayout , unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
void assert_shapes_are_equal(
  const Shape<xLayout,xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
  const Shape<yLayout,ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  typedef Shape<xLayout,xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> x_type ;
  typedef Shape<yLayout,ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> y_type ;

  if ( x != y ) {
    assert_shapes_are_equal_throw(
      typeid(typename x_type::array_layout),
      x_type::value_size ,
      x_type::rank, x.Stride, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7,
      typeid(typename y_type::array_layout),
      y_type::value_size ,
      y_type::rank, y.Stride, y.N0, y.N1, y.N2, y.N3, y.N4, y.N5, y.N6, y.N7 );
  }
}

template< class    xLayout , unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          class    yLayout , unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
void assert_shapes_equal_dimension(
  const Shape<xLayout,xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
  const Shape<yLayout,ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  typedef Shape<xLayout,xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> x_type ;
  typedef Shape<yLayout,ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> y_type ;

  if ( x.rank != y.rank ||
       x.N0 != y.N0 || x.N1 != y.N1 || x.N2 != y.N2 || x.N3 != y.N3 ||
       x.N4 != y.N4 || x.N5 != y.N5 || x.N6 != y.N6 || x.N7 != y.N7 ) {
    assert_shapes_are_equal_throw(
      typeid(typename x_type::array_layout),
      x_type::value_size ,
      x_type::rank, x.Stride, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7,
      typeid(typename y_type::array_layout),
      y_type::value_size ,
      y_type::rank, y.Stride, y.N0, y.N1, y.N2, y.N3, y.N4, y.N5, y.N6, y.N7 );
  }
}

//----------------------------------------------------------------------------

template< class ShapeType > struct assert_shape_is_rank_zero ;
template< class ShapeType > struct assert_shape_is_rank_one ;

template< class Layout , unsigned Size >
struct assert_shape_is_rank_zero< Shape<Layout,Size,0> >
  : public true_type {};

template< class Layout , unsigned Size , unsigned s0 >
struct assert_shape_is_rank_one< Shape<Layout,Size,1,s0> >
  : public true_type {};

//----------------------------------------------------------------------------

void assert_shape_bounds_throw( const size_t rank ,
                                const size_t n0 , const size_t n1 ,
                                const size_t n2 , const size_t n3 ,
                                const size_t n4 , const size_t n5 ,
                                const size_t n6 , const size_t n7 ,
                                const size_t i0 , const size_t i1 ,
                                const size_t i2 , const size_t i3 ,
                                const size_t i4 , const size_t i5 ,
                                const size_t i6 , const size_t i7 );

template< class ShapeType >
inline
void assert_shape_bounds( const ShapeType & shape ,
                          const size_t i0 = 0 ,
                          const size_t i1 = 0 ,
                          const size_t i2 = 0 ,
                          const size_t i3 = 0 ,
                          const size_t i4 = 0 ,
                          const size_t i5 = 0 ,
                          const size_t i6 = 0 ,
                          const size_t i7 = 0 )
{
  const bool ok = 0 == ShapeType::rank ? true : i0 < shape.N0 && (
                  1 == ShapeType::rank ? true : i1 < shape.N1 && (
                  2 == ShapeType::rank ? true : i2 < shape.N2 && (
                  3 == ShapeType::rank ? true : i3 < shape.N3 && (
                  4 == ShapeType::rank ? true : i4 < shape.N4 && (
                  5 == ShapeType::rank ? true : i5 < shape.N5 && (
                  6 == ShapeType::rank ? true : i6 < shape.N6 && (
                  7 == ShapeType::rank ? true : i7 < shape.N7 )))))));

  if ( ! ok ) {
    assert_shape_bounds_throw( ShapeType::rank ,
                               shape.N0 , shape.N1 , shape.N2 , shape.N3 ,
                               shape.N4 , shape.N5 , shape.N6 , shape.N7 ,
                               i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 );
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Specialization and optimization for the Rank 0 shape.

template < class Layout , unsigned ValueSize >
struct Shape< Layout , ValueSize , 0, 1,1,1,1, 1,1,1,1 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 0 ;
  static const unsigned StaticN0     = 1 ;
  static const unsigned StaticN1     = 1 ;
  static const unsigned StaticN2     = 1 ;
  static const unsigned StaticN3     = 1 ;
  static const unsigned StaticN4     = 1 ;
  static const unsigned StaticN5     = 1 ;
  static const unsigned StaticN6     = 1 ;
  static const unsigned StaticN7     = 1 ;

  static const unsigned Stride = 1 ;
  static const unsigned N0 = 1 ;
  static const unsigned N1 = 1 ;
  static const unsigned N2 = 1 ;
  static const unsigned N3 = 1 ;
  static const unsigned N4 = 1 ;
  static const unsigned N5 = 1 ;
  static const unsigned N6 = 1 ;
  static const unsigned N7 = 1 ;

  template< class MemorySpace >
  static inline
  Shape create()
  { return Shape(); }

  static inline
  Shape create( const Shape )
  { return Shape(); }
};

//----------------------------------------------------------------------------
// All-static dimension array

template < class Layout ,
           unsigned ValueSize ,
           unsigned Rank ,
           unsigned s0 , 
           unsigned s1 , 
           unsigned s2 , 
           unsigned s3 , 
           unsigned s4 , 
           unsigned s5 , 
           unsigned s6 , 
           unsigned s7 >
struct Shape {
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = s0 ;
  static const unsigned StaticN1     = s1 ;
  static const unsigned StaticN2     = s2 ;
  static const unsigned StaticN3     = s3 ;
  static const unsigned StaticN4     = s4 ;
  static const unsigned StaticN5     = s5 ;
  static const unsigned StaticN6     = s6 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;

  static const unsigned N0 = s0 ;
  static const unsigned N1 = s1 ;
  static const unsigned N2 = s2 ;
  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create()
  {
    Shape shape ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  static inline
  Shape create( const Shape input )
  {
    Shape shape ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 1 == dynamic_rank <= rank <= 8
template < class Layout ,
           unsigned ValueSize ,
           unsigned Rank ,
           unsigned s1 , 
           unsigned s2 , 
           unsigned s3 , 
           unsigned s4 , 
           unsigned s5 , 
           unsigned s6 , 
           unsigned s7 >
struct Shape< Layout , ValueSize , Rank , 0,s1,s2,s3, s4,s5,s6,s7 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = 0 ;
  static const unsigned StaticN1     = s1 ;
  static const unsigned StaticN2     = s2 ;
  static const unsigned StaticN3     = s3 ;
  static const unsigned StaticN4     = s4 ;
  static const unsigned StaticN5     = s5 ;
  static const unsigned StaticN6     = s6 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;
  unsigned N0 ;

  static const unsigned N1 = s1 ;
  static const unsigned N2 = s2 ;
  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0 , N1, N2, N3, N4, N5, N6, N7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ; // May or may not be static
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 2 == dynamic_rank <= rank <= 8
template < class Layout , unsigned ValueSize , unsigned Rank ,
           unsigned s2 , 
           unsigned s3 , 
           unsigned s4 , 
           unsigned s5 , 
           unsigned s6 , 
           unsigned s7 >
struct Shape< Layout , ValueSize , Rank , 0,0,s2,s3, s4,s5,s6,s7 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 2 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = 0 ;
  static const unsigned StaticN1     = 0 ;
  static const unsigned StaticN2     = s2 ;
  static const unsigned StaticN3     = s3 ;
  static const unsigned StaticN4     = s4 ;
  static const unsigned StaticN5     = s5 ;
  static const unsigned StaticN6     = s6 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;

  static const unsigned N2 = s2 ;
  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 , unsigned SrcN1 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0 , SrcN1, N2, N3, N4, N5, N6, N7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ; // May or may not be static
    shape.N1 = input.N1 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 3 == dynamic_rank <= rank <= 8
template < class Layout , unsigned Rank , unsigned ValueSize ,
           unsigned s3 , 
           unsigned s4 , 
           unsigned s5 , 
           unsigned s6 , 
           unsigned s7 >
struct Shape< Layout , ValueSize , Rank , 0,0,0,s3, s4,s5,s6,s7>
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 3 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = 0 ;
  static const unsigned StaticN1     = 0 ;
  static const unsigned StaticN2     = 0 ;
  static const unsigned StaticN3     = s3 ;
  static const unsigned StaticN4     = s4 ;
  static const unsigned StaticN5     = s5 ;
  static const unsigned StaticN6     = s6 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;

  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 ,
                const unsigned arg_N2 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.N2 = arg_N2 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0,SrcN1,SrcN2,N3,N4,N5,N6,N7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.N2 = input.N2 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 4 == dynamic_rank <= rank <= 8
template < class Layout , unsigned ValueSize , unsigned Rank ,
           unsigned s4 , 
           unsigned s5 , 
           unsigned s6 , 
           unsigned s7 >
struct Shape< Layout , ValueSize , Rank, 0,0,0,0, s4,s5,s6,s7 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 4 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = 0 ;
  static const unsigned StaticN1     = 0 ;
  static const unsigned StaticN2     = 0 ;
  static const unsigned StaticN3     = 0 ;
  static const unsigned StaticN4     = s4 ;
  static const unsigned StaticN5     = s5 ;
  static const unsigned StaticN6     = s6 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;

  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 ,
                const unsigned arg_N2 ,
                const unsigned arg_N3 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.N2 = arg_N2 ;
    shape.N3 = arg_N3 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 , unsigned SrcN3 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0,SrcN1,SrcN2,SrcN3,N4,N5,N6,N7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.N2 = input.N2 ;
    shape.N3 = input.N3 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 5 == dynamic_rank <= rank <= 8
template < class Layout , unsigned ValueSize , unsigned Rank ,
           unsigned s5 , 
           unsigned s6 , 
           unsigned s7 >
struct Shape< Layout , ValueSize , Rank , 0,0,0,0, 0,s5,s6,s7 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 5 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = 0 ;
  static const unsigned StaticN1     = 0 ;
  static const unsigned StaticN2     = 0 ;
  static const unsigned StaticN3     = 0 ;
  static const unsigned StaticN4     = 0 ;
  static const unsigned StaticN5     = s5 ;
  static const unsigned StaticN6     = s6 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;

  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 ,
                const unsigned arg_N2 ,
                const unsigned arg_N3 ,
                const unsigned arg_N4 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.N2 = arg_N2 ;
    shape.N3 = arg_N3 ;
    shape.N4 = arg_N4 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 , unsigned SrcN3 ,
            unsigned SrcN4 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0,SrcN1,SrcN2,SrcN3,SrcN4,N5,N6,N7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.N2 = input.N2 ;
    shape.N3 = input.N3 ;
    shape.N4 = input.N4 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 6 == dynamic_rank <= rank <= 8
template < class Layout , unsigned ValueSize , unsigned Rank ,
           unsigned s6 , 
           unsigned s7 >
struct Shape< Layout , ValueSize , Rank , 0,0,0,0, 0,0,s6,s7 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 6 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = 0 ;
  static const unsigned StaticN1     = 0 ;
  static const unsigned StaticN2     = 0 ;
  static const unsigned StaticN3     = 0 ;
  static const unsigned StaticN4     = 0 ;
  static const unsigned StaticN5     = 0 ;
  static const unsigned StaticN6     = s6 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;

  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 ,
                const unsigned arg_N2 ,
                const unsigned arg_N3 ,
                const unsigned arg_N4 ,
                const unsigned arg_N5 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.N2 = arg_N2 ;
    shape.N3 = arg_N3 ;
    shape.N4 = arg_N4 ;
    shape.N5 = arg_N5 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 , unsigned SrcN3 ,
            unsigned SrcN4 , unsigned SrcN5 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0,SrcN1,SrcN2,SrcN3,SrcN4,SrcN5,N6,N7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.N2 = input.N2 ;
    shape.N3 = input.N3 ;
    shape.N4 = input.N4 ;
    shape.N5 = input.N5 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 7 == dynamic_rank <= rank <= 8
template < class Layout , unsigned ValueSize , unsigned Rank ,
           unsigned s7 >
struct Shape< Layout , ValueSize , Rank , 0,0,0,0, 0,0,0,s7 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 7 ;
  static const unsigned rank         = Rank ;
  static const unsigned StaticN0     = 0 ;
  static const unsigned StaticN1     = 0 ;
  static const unsigned StaticN2     = 0 ;
  static const unsigned StaticN3     = 0 ;
  static const unsigned StaticN4     = 0 ;
  static const unsigned StaticN5     = 0 ;
  static const unsigned StaticN6     = 0 ;
  static const unsigned StaticN7     = s7 ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;

  static const unsigned N7 = s7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 ,
                const unsigned arg_N2 ,
                const unsigned arg_N3 ,
                const unsigned arg_N4 ,
                const unsigned arg_N5 ,
                const unsigned arg_N6 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.N2 = arg_N2 ;
    shape.N3 = arg_N3 ;
    shape.N4 = arg_N4 ;
    shape.N5 = arg_N5 ;
    shape.N6 = arg_N6 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 , unsigned SrcN3 ,
            unsigned SrcN4 , unsigned SrcN5 , unsigned SrcN6 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0,SrcN1,SrcN2,SrcN3,SrcN4,SrcN5,SrcN6,N7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.N2 = input.N2 ;
    shape.N3 = input.N3 ;
    shape.N4 = input.N4 ;
    shape.N5 = input.N5 ;
    shape.N6 = input.N6 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

// 8 == dynamic_rank <= rank <= 8
template < class Layout , unsigned ValueSize >
struct Shape< Layout , ValueSize , 8 , 0,0,0,0, 0,0,0,0 >
{
  typedef Layout  array_layout ;

  static const unsigned value_size   = ValueSize ;
  static const unsigned rank_dynamic = 8 ;
  static const unsigned rank         = 8 ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  unsigned N7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 ,
                const unsigned arg_N2 ,
                const unsigned arg_N3 ,
                const unsigned arg_N4 ,
                const unsigned arg_N5 ,
                const unsigned arg_N6 ,
                const unsigned arg_N7 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.N2 = arg_N2 ;
    shape.N3 = arg_N3 ;
    shape.N4 = arg_N4 ;
    shape.N5 = arg_N5 ;
    shape.N6 = arg_N6 ;
    shape.N7 = arg_N7 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape>::template stride<MemorySpace>( shape );
    return shape ;
  }

  template< unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 , unsigned SrcN3 ,
            unsigned SrcN4 , unsigned SrcN5 , unsigned SrcN6 , unsigned SrcN7 >
  static inline
  Shape create( const Shape< array_layout , value_size , rank ,
                             SrcN0,SrcN1,SrcN2,SrcN3,
                             SrcN4,SrcN5,SrcN6,SrcN7 > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.N2 = input.N2 ;
    shape.N3 = input.N3 ;
    shape.N4 = input.N4 ;
    shape.N5 = input.N5 ;
    shape.N6 = input.N6 ;
    shape.N7 = input.N7 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstShape , class SrcShape >
struct SubShape ;

template< class Layout , unsigned ValueSize , unsigned Rank ,
          unsigned DstN0 , unsigned DstN1 , unsigned DstN2 , unsigned DstN3 ,
          unsigned DstN4 , unsigned DstN5 , unsigned DstN6 , unsigned DstN7 ,
          unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 , unsigned SrcN3 ,
          unsigned SrcN4 , unsigned SrcN5 , unsigned SrcN6 , unsigned SrcN7 >
struct SubShape< Shape< Layout , ValueSize , Rank ,
                        DstN0 , DstN1 , DstN2 , DstN3 ,
                        DstN4 , DstN5 , DstN6 , DstN7 > ,
                 Shape< Layout , ValueSize , Rank ,
                        SrcN0 , SrcN1 , SrcN2 , SrcN3 ,
                        SrcN4 , SrcN5 , SrcN6 , SrcN7 > >
{
  typedef Shape< Layout , ValueSize , Rank ,
                 DstN0 , DstN1 , DstN2 , DstN3 ,
                 DstN4 , DstN5 , DstN6 , DstN7 > DstShape ;

  typedef Shape< Layout , ValueSize , Rank ,
                 SrcN0 , SrcN1 , SrcN2 , SrcN3 ,
                 SrcN4 , SrcN5 , SrcN6 , SrcN7 > SrcShape ;

  typedef typename
    StaticAssert< SrcShape::rank_dynamic <= DstShape::rank_dynamic , SubShape >
      ::type type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src )
  {
    shape  = DstShape::create( src );
    offset = 0 ;
  }
};

//----------------------------------------------------------------------------
// Shape and offset for a single member:

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout , unsigned s0 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 1, s0 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 1, s0 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0);
  }
};

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout , unsigned s0 , unsigned s1 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 2, s0,s1 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 2, s0,s1 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 , const size_t i1 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0,i1);
  }
};

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout , unsigned s0 , unsigned s1 , unsigned s2 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 3, s0,s1,s2 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 3, s0,s1,s2 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 , const size_t i1 ,
                                 const size_t i2 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0,i1,i2);
  }
};

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 4, s0,s1,s2,s3 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 4, s0,s1,s2,s3 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 , const size_t i1 ,
                                 const size_t i2 , const size_t i3 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0,i1,i2,i3);
  }
};

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 5, s0,s1,s2,s3, s4 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 5, s0,s1,s2,s3, s4 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 , const size_t i1 ,
                                 const size_t i2 , const size_t i3 ,
                                 const size_t i4 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0,i1,i2,i3,i4);
  }
};

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 6, s0,s1,s2,s3, s4,s5 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 6, s0,s1,s2,s3, s4,s5 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 , const size_t i1 ,
                                 const size_t i2 , const size_t i3 ,
                                 const size_t i4 , const size_t i5 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0,i1,i2,i3,i4,i5);
  }
};

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 7, s0,s1,s2,s3, s4,s5,s6 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 7, s0,s1,s2,s3, s4,s5,s6 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 , const size_t i1 ,
                                 const size_t i2 , const size_t i3 ,
                                 const size_t i4 , const size_t i5 ,
                                 const size_t i6 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0,i1,i2,i3,i4,i5,i6);
  }
};

template< class DstLayout , unsigned ValueSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 >
struct SubShape< Shape< DstLayout, ValueSize, 0 > ,
                 Shape< SrcLayout, ValueSize, 8, s0,s1,s2,s3, s4,s5,s6,s7 > >
{
  typedef Shape< DstLayout, ValueSize, 0 > DstShape ;
  typedef Shape< SrcLayout, ValueSize, 8, s0,s1,s2,s3, s4,s5,s6,s7 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 , const size_t i1 ,
                                 const size_t i2 , const size_t i3 ,
                                 const size_t i4 , const size_t i5 ,
                                 const size_t i6 , const size_t i7 )
  {
    offset = ShapeMap<SrcShape>::offset(src,i0,i1,i2,i3,i4,i5,i6,i7);
  }
};

//----------------------------------------------------------------------------
// Shape and offset for a span of a rank-1 array.

template< class Layout , unsigned ValueSize >
struct SubShape< Shape< Layout, ValueSize, 1, 0 > ,
                 Shape< Layout, ValueSize, 1, 0 > >
{
  typedef Shape< Layout, ValueSize, 1, 0 > DstShape ;
  typedef Shape< Layout, ValueSize, 1, 0 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src )
  {
    shape  = DstShape::create( src );
    offset = 0 ;
  }

  template< typename iType >
  SubShape( const SrcShape src , const std::pair<iType,iType> span )
  {
    if ( span.first < span.second ) {
      assert_shape_bounds( src , span.first );
      assert_shape_bounds( src , span.second - 1 );
      shape.N0     = span.second - span.first ;
      shape.Stride = src.Stride ;
      offset       = span.first ;
    }
    else {
      shape.N0     = 0 ;
      shape.Stride = 0 ;
      offset       = 0 ;
    }
  }
};

//----------------------------------------------------------------------------

template< class Layout , unsigned ValueSize , unsigned Rank ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 >
size_t cardinality_count(
  const Shape<Layout,ValueSize,Rank,s0,s1,s2,s3,s4,s5,s6,s7> & shape )
{
  return shape.N0 * shape.N1 * shape.N2 * shape.N3 *
         shape.N4 * shape.N5 * shape.N6 * shape.N7 ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_ARRAYSHAPE_HPP */

