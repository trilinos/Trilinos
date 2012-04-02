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

#ifndef KOKKOS_CRSARRAY_HPP
#define KOKKOS_CRSARRAY_HPP

#include <string>
#include <impl/Kokkos_forward.hpp>
#include <impl/Kokkos_ArrayTraits.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Compressed row storage array.
 *
 *  Define a range of entries for each row:
 *    row_entry_begin(row) <= entry < row_entry_end(row)
 *    row_entry_begin(row) == row_entry_end(row-1)
 *
 *  Access data member via 'operator()'
 */
template< class ArrayType , class DeviceType ,
          typename SizeType = typename DeviceType::size_type >
class CrsArray {
public:
  typedef DeviceType  device_type ;
  typedef SizeType    size_type ;
  typedef ArrayType   array_type ;

  typedef typename Impl::remove_all_extents<array_type>::type  value_type ;

  typedef CrsArray< array_type , typename HostMapped< device_type >::type , size_type >
          HostMirror ;

  /** \brief  Rank is ( entry( row , column ) , dimensions<ValueType> ) */
  enum { Rank = Impl::rank<array_type>::value + 1 };

  /*------------------------------------------------------------------*/
  /** \brief  Number of rows */
  size_type row_count() const ;

  /** \brief  Number of entries == row_entry_end( row_count ) */
  size_type entry_count() const ;

  /** \brief  Begining of entries for a given row */
  template< typename iTypeRow >
  size_type row_entry_begin( const iTypeRow & row ) const ;

  /** \brief  End of entries for a given row */
  template< typename iTypeRow >
  size_type row_entry_end(   const iTypeRow & row ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Value for rank-1 array */
  template< typename iTypeEntry >
  inline
  value_type & operator()( const iTypeEntry & entry ) const ;

  /** \brief  Value for rank-2 array */
  template< typename iTypeEntry , typename iType1 >
  inline
  value_type & operator()( const iTypeEntry & entry ,
                           const iType1     & i1 ) const ;

  /** \brief  Value for rank-3 array */
  template< typename iTypeEntry , typename iType1 , typename iType2 >
  inline
  value_type & operator()( const iTypeEntry & entry ,
                           const iType1     & i1 ,
                           const iType2     & i2 ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  CrsArray();

  /** \brief  Construct a view of the array */
  CrsArray( const CrsArray & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  CrsArray & operator = ( const CrsArray & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~CrsArray();

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const CrsArray & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const CrsArray & ) const ;
};

//----------------------------------------------------------------------------

template< class CrsArrayType , class InputType >
inline
CrsArray< typename CrsArrayType::array_type ,
          typename CrsArrayType::device_type ,
          typename CrsArrayType::size_type >
create_crsarray( const std::string & label ,
                 const InputType & input )
{
  return Impl::Factory< CrsArrayType , InputType >::create( label , input );
}

template< class CrsArrayType , class InputType >
inline
CrsArray< typename CrsArrayType::array_type ,
          typename CrsArrayType::device_type ,
          typename CrsArrayType::size_type >
create_crsarray( const InputType & input )
{
  return Impl::Factory< CrsArrayType , InputType >
             ::create( std::string() , input );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/Kokkos_CrsArray_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CRSARRAY_HPP */

