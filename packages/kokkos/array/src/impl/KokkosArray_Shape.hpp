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
#include <impl/KokkosArray_AnalyzeShape.hpp>

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
          class ShapeStaticInfo ,
          unsigned RankDynamic = ShapeStaticInfo::rank_dynamic ,
          unsigned Rank        = ShapeStaticInfo::rank >
struct Shape ;

template< class ShapeType , class MemorySpace >
struct ShapeMap ;

template< class ShapeType >
struct ShapeOffset ;

//----------------------------------------------------------------------------
/** \brief  Shape equality if the value type, layout, and dimensions
 *          are equal.
 */
template <
  class xLayout , class xDataType , unsigned xRankDynamic , unsigned xRank ,
  class yLayout , class yDataType , unsigned yRankDynamic , unsigned yRank >
bool operator == ( const Shape<xLayout,xDataType,xRankDynamic,xRank> & x ,
                   const Shape<yLayout,yDataType,yRankDynamic,yRank> & y )
{
  typedef Shape<xLayout,xDataType,xRankDynamic,xRank> x_type ;
  typedef Shape<yLayout,yDataType,yRankDynamic,yRank> y_type ;

  enum { same_size = x_type::value_size == y_type::value_size };
  enum { same_rank = xRank == yRank };

  // the array layout only matters for 1 < rank
  enum { same_layout = xRank < 2 ||
                       is_same< typename x_type::array_layout ,
                                typename y_type::array_layout >::value };

  return same_size && same_layout && same_rank &&
         x.N0 == y.N0 && x.N1 == y.N1 && x.N2 == y.N2 && x.N3 == y.N3 &&
         x.N4 == y.N4 && x.N5 == y.N5 && x.N6 == y.N6 && x.N7 == y.N7 &&
         x.Stride == y.Stride ;
}

template <
  class xLayout , class xDataType , unsigned xRankDynamic , unsigned xRank ,
  class yLayout , class yDataType , unsigned yRankDynamic , unsigned yRank >
bool operator != ( const Shape<xLayout,xDataType,xRankDynamic,xRank> & x ,
                   const Shape<yLayout,yDataType,yRankDynamic,yRank> & y )
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

template <
  class xLayout , class xDataType , unsigned xRankDynamic , unsigned xRank ,
  class yLayout , class yDataType , unsigned yRankDynamic , unsigned yRank >
inline
void assert_shapes_are_equal(
  const Shape<xLayout,xDataType,xRankDynamic,xRank> & x ,
  const Shape<yLayout,yDataType,yRankDynamic,yRank> & y )
{
  typedef Shape<xLayout,xDataType,xRankDynamic,xRank> x_type ;
  typedef Shape<yLayout,yDataType,yRankDynamic,yRank> y_type ;

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

template <
  class xLayout , class xDataType , unsigned xRankDynamic , unsigned xRank ,
  class yLayout , class yDataType , unsigned yRankDynamic , unsigned yRank >
inline
void assert_shapes_equal_dimension(
  const Shape<xLayout,xDataType,xRankDynamic,xRank> & x ,
  const Shape<yLayout,yDataType,yRankDynamic,yRank> & y )
{
  typedef Shape<xLayout,xDataType,xRankDynamic,xRank> x_type ;
  typedef Shape<yLayout,yDataType,yRankDynamic,yRank> y_type ;

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

void assert_shape_effective_rank1_at_leastN_throw(
  const size_t x_rank , const size_t x_N0 ,
  const size_t x_N1 ,   const size_t x_N2 ,
  const size_t x_N3 ,   const size_t x_N4 ,
  const size_t x_N5 ,   const size_t x_N6 ,
  const size_t x_N7 ,
  const size_t N0 );

template <
  class xLayout , class xDataType , unsigned xRankDynamic , unsigned xRank >
void assert_shape_effective_rank1_at_leastN(
  const Shape<xLayout,xDataType,xRankDynamic,xRank> & x ,
  const size_t n0 )
{
  if ( 1 < x.N1 || 1 < x.N2 || 1 < x.N3 || 1 < x.N4 ||
       1 < x.N5 || 1 < x.N6 || 1 < x.N7 ||
       x.N0 < n0 ) {
     assert_shape_effective_rank1_at_leastN_throw(
       x.rank, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7, n0 );
  }
}

//----------------------------------------------------------------------------

template< class ShapeType > struct assert_shape_is_rank_zero ;
template< class ShapeType > struct assert_shape_is_rank_one ;
template< class ShapeType > struct assert_shape_is_rank_two ;
template< class ShapeType > struct assert_shape_is_rank_three ;
template< class ShapeType > struct assert_shape_is_rank_four ;
template< class ShapeType > struct assert_shape_is_rank_five ;
template< class ShapeType > struct assert_shape_is_rank_six ;
template< class ShapeType > struct assert_shape_is_rank_seven ;
template< class ShapeType > struct assert_shape_is_rank_eight ;

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_zero< Shape< Layout , Type , RankDynamic , 0 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_one< Shape< Layout , Type , RankDynamic , 1 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_two< Shape< Layout , Type , RankDynamic , 2 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_three< Shape< Layout , Type , RankDynamic , 3 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_four< Shape< Layout , Type , RankDynamic , 4 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_five< Shape< Layout , Type , RankDynamic , 5 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_six< Shape< Layout , Type , RankDynamic , 6 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_seven< Shape< Layout , Type , RankDynamic , 7 > >
  : public true_type {};

template < class Layout , class Type , unsigned RankDynamic >
struct assert_shape_is_rank_eight< Shape< Layout , Type , RankDynamic , 8 > >
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

template < class Layout , class Info >
struct Shape< Layout , Info , 0 , 0 >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;
  static const unsigned Stride       = 1 ;

  static const unsigned N0 = Info::N0 ;
  static const unsigned N1 = Info::N1 ;
  static const unsigned N2 = Info::N2 ;
  static const unsigned N3 = Info::N3 ;
  static const unsigned N4 = Info::N4 ;
  static const unsigned N5 = Info::N5 ;
  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

  template< class MemorySpace >
  static inline
  Shape create()
  { return Shape(); }

  static inline
  Shape create( const Shape )
  { return Shape(); }
};

//----------------------------------------------------------------------------

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 0 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;

  static const unsigned N0 = Info::N0 ;
  static const unsigned N1 = Info::N1 ;
  static const unsigned N2 = Info::N2 ;
  static const unsigned N3 = Info::N3 ;
  static const unsigned N4 = Info::N4 ;
  static const unsigned N5 = Info::N5 ;
  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

  template< class MemorySpace >
  static inline
  Shape create()
  {
    Shape shape ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
  {
    Shape shape ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 1 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;
  unsigned N0 ;

  static const unsigned N1 = Info::N1 ;
  static const unsigned N2 = Info::N2 ;
  static const unsigned N3 = Info::N3 ;
  static const unsigned N4 = Info::N4 ;
  static const unsigned N5 = Info::N5 ;
  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 2 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;

  static const unsigned N2 = Info::N2 ;
  static const unsigned N3 = Info::N3 ;
  static const unsigned N4 = Info::N4 ;
  static const unsigned N5 = Info::N5 ;
  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

  template< class MemorySpace >
  static inline
  Shape create( const unsigned arg_N0 ,
                const unsigned arg_N1 )
  {
    Shape shape ;
    shape.N0 = arg_N0 ;
    shape.N1 = arg_N1 ;
    shape.Stride = 0 ; // to suppress compiler warning
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 3 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;

  static const unsigned N3 = Info::N3 ;
  static const unsigned N4 = Info::N4 ;
  static const unsigned N5 = Info::N5 ;
  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

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
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
  {
    Shape shape ;
    shape.N0 = input.N0 ;
    shape.N1 = input.N1 ;
    shape.N2 = input.N2 ;
    shape.Stride = input.Stride ;
    return shape ;
  }
};

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 4 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;

  static const unsigned N4 = Info::N4 ;
  static const unsigned N5 = Info::N5 ;
  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

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
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
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

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 5 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;

  static const unsigned N5 = Info::N5 ;
  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

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
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
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

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 6 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;

  static const unsigned N6 = Info::N6 ;
  static const unsigned N7 = Info::N7 ;

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
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
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

template< class Layout , class Info , unsigned Rank >
struct Shape< Layout , Info , 7 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;

  static const unsigned N7 = Info::N7 ;

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
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  // The source shape may have a lower dynamic rank
  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
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

template< class Layout , class Info >
struct Shape< Layout , Info , 8 , 8 >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = Info::rank_dynamic ;
  static const unsigned rank         = Info::rank ;
  static const unsigned value_size   = Info::value_size ;

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
    shape.Stride = ShapeMap<Shape,MemorySpace>::stride( shape );
    return shape ;
  }

  // The source shape may have a lower dynamic rank
  template< class SrcType , unsigned SrcRankDyn >
  static inline
  Shape create( const Shape< Layout , SrcType , SrcRankDyn , rank > input )
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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class DstShape , class SrcShape >
struct SubShape ;

template< class Layout ,
          class DstDataType , unsigned DstRankDyn ,
          class SrcDataType , unsigned SrcRankDyn , unsigned Rank >
struct SubShape< Shape< Layout , DstDataType , DstRankDyn , Rank > ,
                 Shape< Layout , SrcDataType , SrcRankDyn , Rank > >
{
  typedef Shape< Layout , DstDataType , DstRankDyn , Rank > DstShape ;
  typedef Shape< Layout , SrcDataType , SrcRankDyn , Rank > SrcShape ;

  typedef typename
    StaticAssert< SrcRankDyn <= DstRankDyn , SubShape >::type type ;

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

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 1 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 1 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src , const size_t i0 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0);
  }
};

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 2 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 2 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src ,
            const size_t i0 , const size_t i1 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0,i1);
  }
};

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 3 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 3 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src ,
            const size_t i0 , const size_t i1 ,
            const size_t i2 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0,i1,i2);
  }
};

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 4 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 4 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src ,
            const size_t i0 , const size_t i1 ,
            const size_t i2 , const size_t i3 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0,i1,i2,i3);
  }
};

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 5 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 5 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src ,
            const size_t i0 , const size_t i1 ,
            const size_t i2 , const size_t i3 ,
            const size_t i4 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0,i1,i2,i3,i4);
  }
};

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 6 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 6 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src ,
            const size_t i0 , const size_t i1 ,
            const size_t i2 , const size_t i3 ,
            const size_t i4 , const size_t i5 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0,i1,i2,i3,i4,i5);
  }
};

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 7 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 7 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src ,
            const size_t i0 , const size_t i1 ,
            const size_t i2 , const size_t i3 ,
            const size_t i4 , const size_t i5 ,
            const size_t i6 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0,i1,i2,i3,i4,i5,i6);
  }
};

template< class DstLayout , class DstDataType ,
          class SrcLayout , class SrcDataType , unsigned SrcRankDyn >
struct SubShape< Shape< DstLayout , DstDataType , 0 , 0 > ,
                 Shape< SrcLayout , SrcDataType , SrcRankDyn , 8 > >
{
  typedef Shape< DstLayout , DstDataType , 0 , 0 > DstShape ;
  typedef Shape< SrcLayout , SrcDataType , SrcRankDyn , 8 > SrcShape ;

  typedef SubShape type ;

  DstShape shape ;
  size_t   offset ;

  SubShape( const SrcShape src ,
            const size_t i0 , const size_t i1 ,
            const size_t i2 , const size_t i3 ,
            const size_t i4 , const size_t i5 ,
            const size_t i6 , const size_t i7 )
  {
    offset = ShapeOffset< SrcShape >::apply(src,i0,i1,i2,i3,i4,i5,i6,i7);
  }
};

//----------------------------------------------------------------------------
// Shape and offset for a span of a rank-1 array.

template< class Layout , class DstDataType , class SrcDataType >
struct SubShape< Shape< Layout , DstDataType , 1 , 1 > ,
                 Shape< Layout , SrcDataType , 1 , 1 > >
{
  typedef Shape< Layout , DstDataType , 1 , 1 > DstShape ;
  typedef Shape< Layout , SrcDataType , 1 , 1 > SrcShape ;

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
      offset       = span.first ;
      shape.N0     = span.second - span.first ;
      shape.Stride = src.Stride ;
    }
    else {
      offset       = 0 ;
      shape.N0     = 0 ;
      shape.Stride = 0 ;
    }
  }
};

//----------------------------------------------------------------------------

template< class Layout , class T , unsigned RankDynamic , unsigned Rank >
inline
size_t cardinality_count( const Shape<Layout,T,RankDynamic,Rank> & shape )
{
  return ( 0 == Rank ? 1 : shape.N0 * (
           1 == Rank ? 1 : shape.N1 * (
           2 == Rank ? 1 : shape.N2 * (
           3 == Rank ? 1 : shape.N3 * (
           4 == Rank ? 1 : shape.N4 * (
           5 == Rank ? 1 : shape.N5 * (
           6 == Rank ? 1 : shape.N6 * (
           7 == Rank ? 1 : shape.N7 ))))))));
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_ARRAYSHAPE_HPP */

