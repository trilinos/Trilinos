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

#ifndef KOKKOS_CRSMAP_HPP
#define KOKKOS_CRSMAP_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <iterator>
#include <limits>
#include <impl/Kokkos_forward.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Compressed row storage map.
 *
 *  one-to-one map ( I , J ) -> index
 *  where I in [ 0 .. N ) and J in [ 0 .. M_I )
 *  and 0 <= index < sum( M_I )
 */
template< class DeviceType ,
          typename SizeType = typename DeviceType::size_type >
class CrsMap {
public:
  typedef DeviceType  device_type ;
  typedef SizeType    size_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Domain of first index is [ 0 .. first_count() - 1 ] */
  size_type first_count() const ;

  /** \brief  Domain of second index is
   *          [ 0 .. second_count(first) - 1 ]
   */
  template< typename iType >
  size_type second_count( const iType & first ) const ;

  /** \brief  Total number of values,
   *          only available on the device.
   */
  size_type size() const ;

  /*------------------------------------------------------------------*/
  /** \brief  Map indices */
  template< typename iType , typename jType >
  size_type operator()( const iType & first , const jType & second ) const ;

  typedef std::pair<size_type,size_type> range_type ;

  /** \brief  Range of indices for a given 'first' index.
   *
   *    second_count(first) == range.second - range.first
   *    operator()(first,0) == range.first ;
   */
  template< typename iType >
  range_type range( const iType & first ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  CrsMap();

  /** \brief  Construct a view of the array */
  CrsMap( const CrsMap & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  CrsMap & operator = ( const CrsMap & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~CrsMap();

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const CrsMap & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const CrsMap & ) const ;
};

//----------------------------------------------------------------------------

namespace Impl {

template< class DstType , class SrcType = DstType > class CreateCrsMap ;

} // namespace Impl

//----------------------------------------------------------------------------
/** \brief  Create a CrsMap with the input row sizes.
 *
 *  Constructed array has the following properties:
 *
 *  first_count()  == std::distance( second_size_begin , second_size_end );
 *  second_count(first) == * iter
 *    where  iter = second_size_begin ; std::advance( iter , first );
 */
template< class CrsMapType , typename IteratorType >
typename Impl::CreateCrsMap< CrsMapType >::type
create_labeled_crsmap( const std::string & label ,
                       const IteratorType second_size_begin ,
                       const IteratorType second_size_end )
{
  // Verify the iterator type dereferences to an integer

  typedef std::iterator_traits<IteratorType> traits ;
  typedef typename traits::value_type int_type ;
  enum { iterator_to_integer = std::numeric_limits<int_type>::is_integer };
  enum { OK = Impl::StaticAssert< iterator_to_integer >::value };

  return Impl::CreateCrsMap< CrsMapType >
    ::create( label , second_size_begin , second_size_end );
}

template< class CrsMapType , typename IteratorType >
typename Impl::CreateCrsMap< CrsMapType >::type
create_crsmap( const std::string & label ,
               const IteratorType second_size_begin ,
               const IteratorType second_size_end )
{ return create_labeled_crsmap<CrsMapType>( std::string() , second_size_begin , second_size_end ); }

template< class DstType , class SrcType >
typename Impl::CreateCrsMap< DstType , SrcType >::type
create_crsmap( const SrcType & src )
{ return Impl::CreateCrsMap< DstType , SrcType >::create( src ); }

//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CRSMAP_HPP */

