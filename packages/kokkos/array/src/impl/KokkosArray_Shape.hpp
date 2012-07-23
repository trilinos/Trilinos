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
#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class T , unsigned J = 0 >
struct replace_type_with_size
  : public unsigned_< sizeof(T) > {};

template< class T , unsigned J >
struct replace_type_with_size< const T , J >
  : public unsigned_< sizeof(T) > {};

template< class T >
struct replace_type_with_size< T[] , 0 >
{
  typedef typename replace_type_with_size<T,1>::type type [0] ;
  static const unsigned value = replace_type_with_size<T,1>::value ;
};

template< class T >
struct replace_type_with_size< const T[] , 0 >
{
  typedef typename replace_type_with_size<T,1>::type type [0] ;
  static const unsigned value = replace_type_with_size<T,1>::value ;
};

template< class T , unsigned J >
struct replace_type_with_size< T[0] , J >
{
  typedef typename replace_type_with_size<T,J+1>::type type [0] ;
  static const unsigned value = replace_type_with_size<T,J+1>::value ;
};

template< class T , unsigned J >
struct replace_type_with_size< const T[0] , J >
{
  typedef typename replace_type_with_size<T,J+1>::type type [0] ;
  static const unsigned value = replace_type_with_size<T,J+1>::value ;
};


template< class T , unsigned N , unsigned J >
struct replace_type_with_size< T[N] , J >
{
  typedef typename replace_type_with_size<T,J+1>::type type [N] ;
  static const unsigned value = replace_type_with_size<T,J+1>::value ;
};

template< class T , unsigned N , unsigned J >
struct replace_type_with_size< const T[N] , J >
{
  typedef typename replace_type_with_size<T,J+1>::type type [N] ;
  static const unsigned value = replace_type_with_size<T,J+1>::value ;
};

//----------------------------------------------------------------------------
/** \brief  A Kokkos Array has dynamic dimensions followed by static dimensions.
 *          rank<Type>::value provides the total rank.
 *          rank_dynamic<Type>::value provides the dynamic sub-rank.
 *
 *  Dynamic dimensions are indicated by an empty dimension or dimension of zero.
 *
 *  Example:
 *    Type = double[0][0][27][3]
 *    rank<Type>::value == 4
 *    rank_dynamic<Type>::value == 2
 */
template < class T , unsigned J = 0 >
struct rank_dynamic : public unsigned_< J > {};

template < class T >
struct rank_dynamic< T[] , 0 >
  : public unsigned_< rank_dynamic<T,1>::value > {};

template < class T , unsigned J >
struct rank_dynamic< T[0] , J >
  : public unsigned_< rank_dynamic<T,J+1>::value > {};

template< class T , unsigned N , unsigned J >
struct rank_dynamic< T[N] , J >
  : public unsigned_< J > {};

//----------------------------------------------------------------------------
/** \brief  The shape of a Kokkos Array with dynamic and static dimensions.
 *          Dynamic dimensions are member values and static dimensions are
 *          'static const' values.
 *
 *  The 'Type' must use zero extents instead of empty extents:
 *    Type = change_empty_extent_to_zero_extent<T>::type
 *
 *  The upper bound on the array rank is eight.
 */
template< class Layout , class Type , unsigned RankDynamic , unsigned Rank >
struct Shape ;

template< class Layout , class Type >
struct DefineShape {
  typedef Shape< Layout ,
                 typename replace_type_with_size<Type>::type ,
                 rank_dynamic< Type >::value ,
                 rank<         Type >::value > type ;
};


/** \brief  Map multi-indices.  Callable from host or device code.
 */
template < class ShapeType , class DeviceSpace >
struct ShapeMap ;

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
// Specialization and optimization for the Rank 0 and Rank 1 shapes.

template < class Layout , class Type >
struct Shape< Layout , Type , 0 , 0 >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 0 ;
  static const unsigned value_size   = Type::value ;
  static const unsigned Stride       = 0 ;

  static const unsigned N0 = 0 ;
  static const unsigned N1 = 0 ;
  static const unsigned N2 = 0 ;
  static const unsigned N3 = 0 ;
  static const unsigned N4 = 0 ;
  static const unsigned N5 = 0 ;
  static const unsigned N6 = 0 ;
  static const unsigned N7 = 0 ;
};

template < class Layout , class Type >
struct Shape< Layout , Type , 0 , 1 >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 1 ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;
  static const unsigned Stride       = 0 ;

  static const unsigned N0 = extent<Type,0>::value ;
  static const unsigned N1 = 0 ;
  static const unsigned N2 = 0 ;
  static const unsigned N3 = 0 ;
  static const unsigned N4 = 0 ;
  static const unsigned N5 = 0 ;
  static const unsigned N6 = 0 ;
  static const unsigned N7 = 0 ;
};

template < class Layout , class Type >
struct Shape< Layout , Type , 1 , 1 >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = 1 ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;
  static const unsigned Stride       = 0 ;

  unsigned N0 ;
  static const unsigned N1 = 0 ;
  static const unsigned N2 = 0 ;
  static const unsigned N3 = 0 ;
  static const unsigned N4 = 0 ;
  static const unsigned N5 = 0 ;
  static const unsigned N6 = 0 ;
  static const unsigned N7 = 0 ;
};

//----------------------------------------------------------------------------
// 1 < Rank

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 0 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;

  static const unsigned N0 = extent<Type,0>::value ;
  static const unsigned N1 = extent<Type,1>::value ;
  static const unsigned N2 = extent<Type,2>::value ;
  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 1 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;

  static const unsigned N1 = extent<Type,1>::value ;
  static const unsigned N2 = extent<Type,2>::value ;
  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 2 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 2 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;

  static const unsigned N2 = extent<Type,2>::value ;
  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 3 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 3 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;

  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 4 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 4 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;

  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 5 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 5 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;

  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 6 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 6 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;

  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 7 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 7 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;

  static const unsigned N7 = extent<Type,7>::value ;
};

template< class Layout , class Type , unsigned Rank >
struct Shape< Layout , Type , 8 , Rank >
{
  typedef Layout  array_layout ;

  static const unsigned rank_dynamic = 8 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = remove_all_extents<Type>::type::value ;

  unsigned Stride ;
  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  unsigned N7 ;
};

//----------------------------------------------------------------------------
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

template< class T , unsigned RankDynamic , unsigned Rank >
inline
size_t allocation_count( const Shape<LayoutLeft,T,RankDynamic,Rank> & shape )
{
  return ( 0 == Rank ? 1 : (
           1 == Rank ? shape.N0 : shape.Stride * shape.N1 * (
           2 == Rank ? 1 : shape.N2 * (
           3 == Rank ? 1 : shape.N3 * (
           4 == Rank ? 1 : shape.N4 * (
           5 == Rank ? 1 : shape.N5 * (
           6 == Rank ? 1 : shape.N6 * (
           7 == Rank ? 1 : shape.N7 ))))))));
}

//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,0,0>, MemorySpace> {

  typedef Shape<LayoutLeft,T,0,0> output_type ;

  inline static
  output_type create()
  {
    output_type shape ;
    return shape ;
  }
};

template < class T , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,0,1>, MemorySpace> {

  typedef Shape<LayoutLeft,T,0,1> output_type ;

  inline static
  output_type create()
  {
    output_type shape ;
    return shape ;
  }
};

template < class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,0,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,0,Rank> output_type ;

  inline static
  output_type create()
  {
    output_type shape ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,1,1>, MemorySpace> {

  typedef Shape<LayoutLeft,T,1,1> output_type ;

  inline static
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    return shape ;
  }
};

template < class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,1,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,1,Rank> output_type ;

  inline static
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

//----------------------------------------------------------------------------

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,2,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,2,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,3,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,3,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,4,Rank>, MemorySpace> {

  typedef Shape<LayoutLeft,T,4,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,5,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,5,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,6,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,6,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,7,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,7,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutLeft,T,8,Rank> , MemorySpace> {

  typedef Shape<LayoutLeft,T,8,Rank> output_type ;

  inline static
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;
    shape.N7 = n7 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ?  MemorySpace::preferred_alignment( shape.value_size , shape.N0 )
      :  shape.N0 ;

    return shape ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class T , unsigned RankDynamic , unsigned Rank >
inline
size_t allocation_count( const Shape<LayoutRight,T,RankDynamic,Rank> & shape )
{
  return ( 0 == Rank ? 1 : shape.N0 * (
           1 == Rank ? 1 : shape.Stride ));
}


template< class ShapeType >
inline
size_t shape_right_block_count( const ShapeType & shape )
{
  enum { Rank = ShapeType::rank };

  return shape.N1 * ( 2 == Rank ? 1 : shape.N2 * (
                      3 == Rank ? 1 : shape.N3 * (
                      4 == Rank ? 1 : shape.N4 * (
                      5 == Rank ? 1 : shape.N5 * (
                      6 == Rank ? 1 : shape.N6 * (
                      7 == Rank ? 1 : shape.N7 ))))));
}

//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory<Shape<LayoutRight,T,0,0>,MemorySpace> {

  typedef Shape<LayoutRight,T,0,0> output_type ;

  inline static
  output_type create()
  {
    return output_type();
  }
};

template < class T , class MemorySpace >
struct Factory<Shape<LayoutRight,T,0,1>,MemorySpace> {

  typedef Shape<LayoutRight,T,0,1> output_type ;

  inline static
  output_type create()
  {
    return output_type();
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,0,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,0,Rank> output_type ;

  inline static 
  output_type create()
  {
    output_type shape ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

//----------------------------------------------------------------------------

template < class T , class MemorySpace >
struct Factory<Shape<LayoutRight,T,1,1>,MemorySpace> {

  typedef Shape<LayoutRight,T,1,1> output_type ;

  inline static 
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,1,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,1,Rank> output_type ;

  inline static 
  output_type create( size_t n0 )
  {
    output_type shape ;
    shape.N0 = n0 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,2,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,2,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,3,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,3,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,4,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,4,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,5,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,5,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,6,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,6,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,7,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,7,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

template< class T , unsigned Rank , class MemorySpace >
struct Factory< Shape<LayoutRight,T,8,Rank> , MemorySpace > {

  typedef Shape<LayoutRight,T,8,Rank> output_type ;

  inline static 
  output_type create( size_t n0 , size_t n1 , size_t n2 , size_t n3 ,
                      size_t n4 , size_t n5 , size_t n6 , size_t n7 )
  {
    output_type shape ;
    shape.N0 = n0 ;
    shape.N1 = n1 ;
    shape.N2 = n2 ;
    shape.N3 = n3 ;
    shape.N4 = n4 ;
    shape.N5 = n5 ;
    shape.N6 = n6 ;
    shape.N7 = n7 ;

    shape.Stride =
      ( 1 < shape.N1 || 1 < shape.N2 || 1 < shape.N3 || 1 < shape.N4 ||
        1 < shape.N5 || 1 < shape.N6 || 1 < shape.N7 )
      ? MemorySpace::preferred_alignment(
          shape.value_size , shape_right_block_count( shape ) )
        : 1 ;

    return shape ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_ARRAYSHAPE_HPP */

