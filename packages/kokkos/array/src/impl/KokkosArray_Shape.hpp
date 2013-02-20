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

#ifndef KOKKOSARRAY_SHAPE_HPP
#define KOKKOSARRAY_SHAPE_HPP

#include <typeinfo>
#include <utility>
#include <KokkosArray_Macros.hpp>
#include <KokkosArray_Layout.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  The shape of a KokkosArray with dynamic and static dimensions.
 *          Dynamic dimensions are member values and static dimensions are
 *          'static const' values.
 *
 *  The upper bound on the array rank is eight.
 */
template< unsigned ScalarSize ,
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

template< class ShapeType , class Layout >
struct ShapeMap ;

//----------------------------------------------------------------------------
/** \brief  Shape equality if the value type, layout, and dimensions
 *          are equal.
 */
template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
bool operator == ( const Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
                   const Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  enum { same_size = xSize == ySize };
  enum { same_rank = xRank == yRank };

  return same_size && same_rank &&
         x.N0 == y.N0 && x.N1 == y.N1 && x.N2 == y.N2 && x.N3 == y.N3 &&
         x.N4 == y.N4 && x.N5 == y.N5 && x.N6 == y.N6 && x.N7 == y.N7 ;
}

template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize ,unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
bool operator != ( const Shape<xSize,xRank,
                               xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
                   const Shape<ySize,yRank,
                               yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{ return ! operator == ( x , y ); }

//----------------------------------------------------------------------------

void assert_shapes_are_equal_throw(
  const unsigned x_scalar_size ,
  const unsigned x_rank ,
  const unsigned x_N0 , const unsigned x_N1 ,
  const unsigned x_N2 , const unsigned x_N3 ,
  const unsigned x_N4 , const unsigned x_N5 ,
  const unsigned x_N6 , const unsigned x_N7 ,

  const unsigned y_scalar_size ,
  const unsigned y_rank ,
  const unsigned y_N0 , const unsigned y_N1 ,
  const unsigned y_N2 , const unsigned y_N3 ,
  const unsigned y_N4 , const unsigned y_N5 ,
  const unsigned y_N6 , const unsigned y_N7 );

template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
void assert_shapes_are_equal(
  const Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
  const Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  typedef Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> x_type ;
  typedef Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> y_type ;

  if ( x != y ) {
    assert_shapes_are_equal_throw(
      x_type::scalar_size, x_type::rank, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7,
      y_type::scalar_size, y_type::rank, y.N0, y.N1, y.N2, y.N3, y.N4, y.N5, y.N6, y.N7 );
  }
}

template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
void assert_shapes_equal_dimension(
  const Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
  const Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  typedef Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> x_type ;
  typedef Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> y_type ;

  if ( x.rank != y.rank ||
       x.N0 != y.N0 || x.N1 != y.N1 || x.N2 != y.N2 || x.N3 != y.N3 ||
       x.N4 != y.N4 || x.N5 != y.N5 || x.N6 != y.N6 || x.N7 != y.N7 ) {
    assert_shapes_are_equal_throw(
      x_type::scalar_size ,
      x_type::rank, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7,
      y_type::scalar_size ,
      y_type::rank, y.N0, y.N1, y.N2, y.N3, y.N4, y.N5, y.N6, y.N7 );
  }
}

//----------------------------------------------------------------------------

template< class ShapeType > struct assert_shape_is_rank_zero ;
template< class ShapeType > struct assert_shape_is_rank_one ;

template< unsigned Size >
struct assert_shape_is_rank_zero< Shape<Size,0> >
  : public true_type {};

template< unsigned Size , unsigned s0 >
struct assert_shape_is_rank_one< Shape<Size,1,s0> >
  : public true_type {};

//----------------------------------------------------------------------------

/** \brief  Array bounds assertion templated on the execution space
 *          to allow device-specific abort code.
 */
template< class ExecutionSpace >
struct AssertShapeBoundsAbort ;

template<>
struct AssertShapeBoundsAbort< KokkosArray::HostSpace >
{
  static void apply( const size_t rank ,
                     const size_t n0 , const size_t n1 ,
                     const size_t n2 , const size_t n3 ,
                     const size_t n4 , const size_t n5 ,
                     const size_t n6 , const size_t n7 ,
                     const size_t i0 , const size_t i1 ,
                     const size_t i2 , const size_t i3 ,
                     const size_t i4 , const size_t i5 ,
                     const size_t i6 , const size_t i7 );
};

template< class ExecutionDevice >
struct AssertShapeBoundsAbort
{
  KOKKOSARRAY_INLINE_FUNCTION
  static void apply( const size_t rank ,
                     const size_t n0 , const size_t n1 ,
                     const size_t n2 , const size_t n3 ,
                     const size_t n4 , const size_t n5 ,
                     const size_t n6 , const size_t n7 ,
                     const size_t i0 , const size_t i1 ,
                     const size_t i2 , const size_t i3 ,
                     const size_t i4 , const size_t i5 ,
                     const size_t i6 , const size_t i7 )
    {
      AssertShapeBoundsAbort< KokkosArray::HostSpace >
        ::apply( rank , n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 ,
                        i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 );
    }
};

template< class ShapeType >
KOKKOSARRAY_INLINE_FUNCTION
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
    AssertShapeBoundsAbort< ExecutionSpace >
      ::apply( ShapeType::rank ,
               shape.N0 , shape.N1 , shape.N2 , shape.N3 ,
               shape.N4 , shape.N5 , shape.N6 , shape.N7 ,
               i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 );
  }
}

#if defined( KOKKOSARRAY_EXPRESSION_CHECK )
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( S , I0 ) assert_shape_bounds(S,I0);
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( S , I0 , I1 ) assert_shape_bounds(S,I0,I1);
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( S , I0 , I1 , I2 ) assert_shape_bounds(S,I0,I1,I2);
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( S , I0 , I1 , I2 , I3 ) assert_shape_bounds(S,I0,I1,I2,I3);
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( S , I0 , I1 , I2 , I3 , I4 ) assert_shape_bounds(S,I0,I1,I2,I3,I4);
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( S , I0 , I1 , I2 , I3 , I4 , I5 ) assert_shape_bounds(S,I0,I1,I2,I3,I4,I5);
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 ) assert_shape_bounds(S,I0,I1,I2,I3,I4,I5,I6);
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 , I7 ) assert_shape_bounds(S,I0,I1,I2,I3,I4,I5,I6,I7);
#else
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( S , I0 ) /* */
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( S , I0 , I1 ) /* */
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( S , I0 , I1 , I2 ) /* */
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( S , I0 , I1 , I2 , I3 ) /* */
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( S , I0 , I1 , I2 , I3 , I4 ) /* */
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( S , I0 , I1 , I2 , I3 , I4 , I5 ) /* */
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 ) /* */
#define KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 , I7 ) /* */
#endif


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Specialization and optimization for the Rank 0 shape.

template < unsigned ScalarSize >
struct Shape< ScalarSize , 0, 1,1,1,1, 1,1,1,1 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 0 ;

  static const unsigned N0 = 1 ;
  static const unsigned N1 = 1 ;
  static const unsigned N2 = 1 ;
  static const unsigned N3 = 1 ;
  static const unsigned N4 = 1 ;
  static const unsigned N5 = 1 ;
  static const unsigned N6 = 1 ;
  static const unsigned N7 = 1 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  {}
};

//----------------------------------------------------------------------------
// All-static dimension array

template < unsigned ScalarSize ,
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

  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = Rank ;

  static const unsigned N0 = s0 ;
  static const unsigned N1 = s1 ;
  static const unsigned N2 = s2 ;
  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  {}
};

// 1 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize ,
           unsigned Rank ,
           unsigned s1 ,
           unsigned s2 ,
           unsigned s3 ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,s1,s2,s3, s4,s5,s6,s7 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = Rank ;

  unsigned N0 ;

  static const unsigned N1 = s1 ;
  static const unsigned N2 = s2 ;
  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned = 0 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; }
};

// 2 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s2 ,
           unsigned s3 ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,s2,s3, s4,s5,s6,s7 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 2 ;
  static const unsigned rank         = Rank ;

  unsigned N0 ;
  unsigned N1 ;

  static const unsigned N2 = s2 ;
  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; }
};

// 3 == dynamic_rank <= rank <= 8
template < unsigned Rank , unsigned ScalarSize ,
           unsigned s3 ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,s3, s4,s5,s6,s7>
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 3 ;
  static const unsigned rank         = Rank ;

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;

  static const unsigned N3 = s3 ;
  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; }
};

// 4 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank, 0,0,0,0, s4,s5,s6,s7 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 4 ;
  static const unsigned rank         = Rank ;

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;

  static const unsigned N4 = s4 ;
  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ; }
};

// 5 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,0, 0,s5,s6,s7 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 5 ;
  static const unsigned rank         = Rank ;

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;

  static const unsigned N5 = s5 ;
  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ; s.N4 = n4 ; }
};

// 6 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,0, 0,0,s6,s7 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 6 ;
  static const unsigned rank         = Rank ;

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;

  static const unsigned N6 = s6 ;
  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned n5 = 0 , unsigned = 0 , unsigned = 0 )
  {
    s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ;
    s.N4 = n4 ; s.N5 = n5 ;
  }
};

// 7 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,0, 0,0,0,s7 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 7 ;
  static const unsigned rank         = Rank ;

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;

  static const unsigned N7 = s7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned n5 , unsigned n6 , unsigned = 0 )
  {
    s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ;
    s.N4 = n4 ; s.N5 = n5 ; s.N6 = n6 ;
  }
};

// 8 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize >
struct Shape< ScalarSize , 8 , 0,0,0,0, 0,0,0,0 >
{
  static const unsigned scalar_size   = ScalarSize ;
  static const unsigned rank_dynamic = 8 ;
  static const unsigned rank         = 8 ;

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  unsigned N7 ;

  KOKKOSARRAY_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7 )
  {
    s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ;
    s.N4 = n4 ; s.N5 = n5 ; s.N6 = n6 ; s.N7 = n7 ;
  }
};

//----------------------------------------------------------------------------

template< class ShapeType , unsigned N ,
          unsigned R = ShapeType::rank_dynamic >
struct ShapeInsert ;

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 0 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 N ,
                 ShapeType::N0 ,
                 ShapeType::N1 ,
                 ShapeType::N2 ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 1 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 N ,
                 ShapeType::N1 ,
                 ShapeType::N2 ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 2 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N2 ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 3 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 4 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 5 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 6 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 7 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N > type ;
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstLayout , class DstShape ,
          class SrcLayout , class SrcShape >
struct SubShape ;

// Equal rank:

template< class Layout ,
          unsigned ScalarSize , unsigned Rank ,
          unsigned DstN0 , unsigned DstN1 , unsigned DstN2 , unsigned DstN3 ,
          unsigned DstN4 , unsigned DstN5 , unsigned DstN6 , unsigned DstN7 ,
          unsigned SrcN0 , unsigned SrcN1 , unsigned SrcN2 , unsigned SrcN3 ,
          unsigned SrcN4 , unsigned SrcN5 , unsigned SrcN6 , unsigned SrcN7 >
struct SubShape< Layout ,
                 Shape< ScalarSize , Rank ,
                        DstN0 , DstN1 , DstN2 , DstN3 ,
                        DstN4 , DstN5 , DstN6 , DstN7 > ,
                 Layout ,
                 Shape< ScalarSize , Rank ,
                        SrcN0 , SrcN1 , SrcN2 , SrcN3 ,
                        SrcN4 , SrcN5 , SrcN6 , SrcN7 > >
{
  typedef Shape< ScalarSize , Rank ,
                 DstN0 , DstN1 , DstN2 , DstN3 ,
                 DstN4 , DstN5 , DstN6 , DstN7 > DstShape ;

  typedef Shape< ScalarSize , Rank ,
                 SrcN0 , SrcN1 , SrcN2 , SrcN3 ,
                 SrcN4 , SrcN5 , SrcN6 , SrcN7 > SrcShape ;

  typedef typename
    StaticAssert< SrcShape::rank_dynamic <= DstShape::rank_dynamic , SubShape >
      ::type type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride )
  {
    DstShape::assign( shape , src.N0 , src.N1 , src.N2 , src.N3 , src.N4 , src.N5 , src.N6 , src.N7 );
    stride = src_stride ;
    offset = 0 ;
  }
};

//----------------------------------------------------------------------------
// Shape and offset for a single member:

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout , unsigned s0 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 1, s0 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 1, s0 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride , const size_t i0 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0);
  }
};

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout , unsigned s0 , unsigned s1 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 2, s0,s1 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 2, s0,s1 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride , const size_t i0 , const size_t i1 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0,i1);
  }
};

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout , unsigned s0 , unsigned s1 , unsigned s2 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 3, s0,s1,s2 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 3, s0,s1,s2 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride ,
            const size_t i0 , const size_t i1 , const size_t i2 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0,i1,i2);
  }
};

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 4, s0,s1,s2,s3 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 4, s0,s1,s2,s3 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride ,
            const size_t i0 , const size_t i1 , const size_t i2 , const size_t i3 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0,i1,i2,i3);
  }
};

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 5, s0,s1,s2,s3, s4 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 5, s0,s1,s2,s3, s4 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride ,
            const size_t i0 , const size_t i1 , const size_t i2 ,
            const size_t i3 , const size_t i4 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0,i1,i2,i3,i4);
  }
};

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 6, s0,s1,s2,s3, s4,s5 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 6, s0,s1,s2,s3, s4,s5 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride ,
            const size_t i0 , const size_t i1 , const size_t i2 , const size_t i3 ,
            const size_t i4 , const size_t i5 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0,i1,i2,i3,i4,i5);
  }
};

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 7, s0,s1,s2,s3, s4,s5,s6 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 7, s0,s1,s2,s3, s4,s5,s6 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride ,
            const size_t i0 , const size_t i1 , const size_t i2 , const size_t i3 ,
            const size_t i4 , const size_t i5 , const size_t i6 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0,i1,i2,i3,i4,i5,i6);
  }
};

template< class DstLayout , unsigned ScalarSize ,
          class SrcLayout ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 >
struct SubShape< DstLayout , Shape< ScalarSize, 0 > ,
                 SrcLayout , Shape< ScalarSize, 8, s0,s1,s2,s3, s4,s5,s6,s7 > >
{
  typedef Shape< ScalarSize, 0 > DstShape ;
  typedef Shape< ScalarSize, 8, s0,s1,s2,s3, s4,s5,s6,s7 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride ,
            const size_t i0 , const size_t i1 , const size_t i2 , const size_t i3 ,
            const size_t i4 , const size_t i5 , const size_t i6 , const size_t i7 )
  {
    stride = 0 ;
    offset = ShapeMap<SrcShape,SrcLayout>::offset(src,src_stride,i0,i1,i2,i3,i4,i5,i6,i7);
  }
};

//----------------------------------------------------------------------------
// Shape and offset for a span of a rank-1 array.

template< class Layout , unsigned ScalarSize >
struct SubShape< Layout , Shape< ScalarSize, 1, 0 > ,
                 Layout , Shape< ScalarSize, 1, 0 > >
{
  typedef Shape< ScalarSize, 1, 0 > DstShape ;
  typedef Shape< ScalarSize, 1, 0 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   stride ;
  size_t   offset ;

  SubShape( const SrcShape src , const unsigned src_stride )
  {
    shape  = src ;
    stride = 1 ;
    offset = 0 ;
  }

  template< typename iType >
  SubShape( const SrcShape src , const unsigned src_stride , const std::pair<iType,iType> span )
  {
    if ( span.first < span.second ) {
      assert_shape_bounds( src , span.first );
      assert_shape_bounds( src , span.second - 1 );
      shape.N0     = span.second - span.first ;
      stride       = 1 ;
      offset       = span.first ;
    }
    else {
      shape.N0     = 0 ;
      stride       = 1 ;
      offset       = 0 ;
    }
  }
};

//----------------------------------------------------------------------------

template< unsigned ScalarSize , unsigned Rank ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 >
size_t cardinality_count(
  const Shape<ScalarSize,Rank,s0,s1,s2,s3,s4,s5,s6,s7> & shape )
{
  return shape.N0 * shape.N1 * shape.N2 * shape.N3 *
         shape.N4 * shape.N5 * shape.N6 * shape.N7 ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_ARRAYSHAPE_HPP */

