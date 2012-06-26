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
template < class X ,
           class T = typename change_empty_extent_to_zero_extent<X>::type >
struct rank_dynamic : public unsigned_<0> {};

template < class X , class T >
struct rank_dynamic< X , T[0] >
  : public unsigned_< rank_dynamic<T>::value + 1 > {};

template< class X , class T , unsigned N >
struct rank_dynamic< X , T[N] >
  : public unsigned_< rank_dynamic<T>::value > {};

//----------------------------------------------------------------------------
/** \brief  Remove the Ith extent */
template < class X , unsigned I ,
           class T = typename change_empty_extent_to_zero_extent<X>::type >
struct remove_extent {
  typedef X type ;
};

template < class X , class T , unsigned N >
struct remove_extent<X,0,T[N]> {
  typedef T type ;
};

template < class X , class T >
struct remove_extent<X,0,T[0]> {
  typedef T type ;
};

template < class X , unsigned I , class T , unsigned N >
struct remove_extent<X,I,T[N]> {
  typedef typename remove_extent<T,I-1>::type type[N] ;
};

template < class X , unsigned I , class T >
struct remove_extent<X,I,T[0]> {
  typedef typename remove_extent<T,I-1>::type type[0] ;
};

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
template< class Layout ,
          class Type ,
          unsigned RankDynamic = rank_dynamic<Type>::value ,
          unsigned Rank        = rank<Type>::value >
struct Shape ;

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

  enum { same_type = is_same< typename x_type::value_type ,
                              typename y_type::value_type >::value };

  enum { same_rank = xRank == yRank };

  // the array layout only matters for 1 < rank
  enum { same_layout = xRank < 2 ||
                       is_same< typename x_type::array_layout ,
                                typename y_type::array_layout >::value };

  return same_type && same_layout && same_rank &&
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
  const std::type_info & x_value_type ,
  const unsigned x_rank , const unsigned x_stride ,
  const unsigned x_N0 , const unsigned x_N1 ,
  const unsigned x_N2 , const unsigned x_N3 ,
  const unsigned x_N4 , const unsigned x_N5 ,
  const unsigned x_N6 , const unsigned x_N7 ,

  const std::type_info & y_layout ,
  const std::type_info & y_value_type ,
  const unsigned y_rank , const unsigned y_stride ,
  const unsigned y_N0 , const unsigned y_N1 ,
  const unsigned y_N2 , const unsigned y_N3 ,
  const unsigned y_N4 , const unsigned y_N5 ,
  const unsigned y_N6 , const unsigned y_N7 );

template <
  class xLayout , class xDataType , unsigned xRankDynamic , unsigned xRank ,
  class yLayout , class yDataType , unsigned yRankDynamic , unsigned yRank >
void assert_shapes_are_equal(
  const Shape<xLayout,xDataType,xRankDynamic,xRank> & x ,
  const Shape<yLayout,yDataType,yRankDynamic,yRank> & y )
{
  typedef Shape<xLayout,xDataType,xRankDynamic,xRank> x_type ;
  typedef Shape<yLayout,yDataType,yRankDynamic,yRank> y_type ;

  if ( x != y ) {
    assert_shapes_are_equal_throw(
      typeid(typename x_type::array_layout),
      typeid(typename x_type::value_type),
      x_type::rank, x.Stride, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7,
      typeid(typename y_type::array_layout),
      typeid(typename y_type::value_type),
      y_type::rank, y.Stride, y.N0, y.N1, y.N2, y.N3, y.N4, y.N5, y.N6, y.N7 );
  }
}

//----------------------------------------------------------------------------

template < class ShapeType , unsigned Rank >
struct assert_shape_is_rank ;

template < class Layout , class Type , unsigned RankDynamic , unsigned Rank >
struct assert_shape_is_rank< Shape< Layout , Type , RankDynamic , Rank > , Rank >
  : public true_type {};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Specialization and optimization for the Rank 0 and Rank 1 shapes.

template < class Layout , class Type >
struct Shape< Layout , Type , 0 , 0 >
{
  typedef Layout array_layout ;
  typedef Type   data_type ;
  typedef Type   value_type ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 0 ;
  static const unsigned value_size   = sizeof(value_type);
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

template < class Layout , class Type , unsigned N >
struct Shape< Layout , Type[N] , 0 , 1 >
{
  typedef Layout  array_layout ;
  typedef Type    data_type[N] ;
  typedef Type    value_type ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 1 ;
  static const unsigned value_size   = sizeof(value_type);
  static const unsigned Stride       = 0 ;

  static const unsigned N0 = N ;
  static const unsigned N1 = 0 ;
  static const unsigned N2 = 0 ;
  static const unsigned N3 = 0 ;
  static const unsigned N4 = 0 ;
  static const unsigned N5 = 0 ;
  static const unsigned N6 = 0 ;
  static const unsigned N7 = 0 ;
};

template < class Layout , class Type >
struct Shape< Layout , Type[0] , 1 , 1 >
{
  typedef Layout  array_layout ;
  typedef Type    data_type[0] ;
  typedef Type    value_type ;

  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = 1 ;
  static const unsigned value_size   = sizeof(value_type);
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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 2 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 3 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 4 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 5 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 6 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 7 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  typedef Layout                                   array_layout ;
  typedef Type                                     data_type ;
  typedef typename remove_all_extents<Type>::type  value_type ;

  static const unsigned rank_dynamic = 8 ;
  static const unsigned rank         = Rank ;
  static const unsigned value_size   = sizeof(value_type);

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
  return ( 0 == Rank ? 1 : shape.Stride * (
           1 == Rank ? 1 : shape.N1 * (
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , shape.N0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , shape.N0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , n0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , n0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , n0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , n0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , n0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , n0 );
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
    shape.Stride = MemorySpace::preferred_alignment( shape.value_size , n0 );
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
    output_type shape ;
    return shape ;
  }
};

template < class T , class MemorySpace >
struct Factory<Shape<LayoutRight,T,0,1>,MemorySpace> {

  typedef Shape<LayoutRight,T,0,1> output_type ;

  inline static
  output_type create()
  {
    output_type shape ;
    return shape ;
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
      MemorySpace::preferred_alignment(
        shape.value_size ,
        shape.N1 * ( 2 == Rank ? 1 : shape.N2 * (
                     3 == Rank ? 1 : shape.N3 * (
                     4 == Rank ? 1 : shape.N4 * (
                     5 == Rank ? 1 : shape.N5 * (
                     6 == Rank ? 1 : shape.N6 * (
                     7 == Rank ? 1 : shape.N7 )))))) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
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
      MemorySpace::preferred_alignment(
        shape.value_size , shape_right_block_count( shape ) );
    return shape ;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_ARRAYSHAPE_HPP */

