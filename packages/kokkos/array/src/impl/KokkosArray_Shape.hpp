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

#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>

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
template < class Type >
struct rank_dynamic ;

/** \brief  The shape of a Kokkos Array with dynamic and static dimensions.
 *          Dynamic dimensions are member values and static dimensions are
 *          'static const' values.
 *
 *  The upper bound on the array rank is eight.
 */
template< class Type ,
          class Layout ,
          unsigned RankDynamic  = rank_dynamic<Type>::value ,
          unsigned RankComplete = rank<Type>::value >
struct Shape ;

/** \brief  Assert that the rank of the array shape is the given value */
template< class ShapeType , unsigned Rank >
struct ShapeAssertRank
  : public StaticAssert< ShapeType::rank == Rank , true_type >::type {};

/** \brief  Query a dimension of an array from its shape.
 *          The dimension may be dynamic or static.
 */
template < class Type , class Layout >
unsigned dimension( const Shape<Type,Layout> & shape , const unsigned i );

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template < class T > struct rank_dynamic      : public unsigned_<0> {};
template < class T > struct rank_dynamic<T[]> : public unsigned_<1> {};
template < class T > struct rank_dynamic<T[0]>
  : public unsigned_< 1 + rank_dynamic<T>::value > {};

//----------------------------------------------------------------------------
// Specialization and optimization for the Rank 0 and Rank 1 shapes.

template < class Type , class Layout >
struct Shape< Type , Layout , 0 , 0 >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 0 ;
  static const unsigned value_size   = sizeof(value_type);
};

template < class Type , class Layout >
struct Shape< Type , Layout , 0 , 1 >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = 1 ;
  static const unsigned value_size   = sizeof(value_type);

  static const unsigned N0 = extent<Type,0>::value ;
};

template < class Type , class Layout >
struct Shape< Type , Layout , 1 , 1 >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = 1 ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
};

//----------------------------------------------------------------------------

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 0 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 0 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  static const unsigned N0 = extent<Type,0>::value ;
  static const unsigned N1 = extent<Type,1>::value ;
  static const unsigned N2 = extent<Type,2>::value ;
  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 1 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 1 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  static const unsigned N1 = extent<Type,1>::value ;
  static const unsigned N2 = extent<Type,2>::value ;
  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 2 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 2 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  unsigned N1 ;
  static const unsigned N2 = extent<Type,2>::value ;
  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 3 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 3 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  static const unsigned N3 = extent<Type,3>::value ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 4 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 4 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  static const unsigned N4 = extent<Type,4>::value ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 5 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 5 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  static const unsigned N5 = extent<Type,5>::value ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 6 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 6 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  static const unsigned N6 = extent<Type,6>::value ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 7 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 7 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  static const unsigned N7 = extent<Type,7>::value ;

  unsigned Stride ;
};

template< class Type , class Layout , unsigned RankComplete >
struct Shape< Type , Layout , 8 , RankComplete >
{
  typedef typename remove_all_extents<Type>::type value_type ;
  typedef Layout                                  layout ;

  static const unsigned rank_dynamic = 8 ;
  static const unsigned rank         = RankComplete ;
  static const unsigned value_size   = sizeof(value_type);

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  unsigned N7 ;

  unsigned Stride ;
};

//----------------------------------------------------------------------------

template < class Type , class Layout >
inline
unsigned dimension( const Shape<Type,Layout> & shape , const unsigned i )
{
  return 0u == i ? shape.N0 : (
         1u == i ? shape.N1 : (
         2u == i ? shape.N2 : (
         3u == i ? shape.N3 : (
         4u == i ? shape.N4 : (
         5u == i ? shape.N5 : (
         6u == i ? shape.N6 : (
         7u == i ? shape.N7 : 0 )))))));
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_ARRAYSHAPE_HPP */

