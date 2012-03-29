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

#include <string>
#include <impl/Kokkos_forward.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

namespace Kokkos {

template< class Device , typename SizeType > class CrsColumnMap ;
template< class Device , typename SizeType > class CrsColumnIdentity ;

//----------------------------------------------------------------------------
/** \brief  Compressed row storage map.
 *
 *  Define a range of entries for each row:
 *    row_entry_begin(row) <= entry < row_entry_end(row)
 *    row_entry_end(row-1) == row_entry_begin(row)
 *
 *  Optionally include a map of entry -> column defining the existence of
 *    ( row , column(entry) )
 */
template< class DeviceType ,
          template< class , typename > class Column = CrsColumnIdentity ,
          typename SizeType = typename DeviceType::size_type >
class CrsMap {
public:
  typedef DeviceType  device_type ;
  typedef SizeType    size_type ;
  typedef CrsMap< size_type , Column , void /* Host */ > HostMirror ;

  /*------------------------------------------------------------------*/
  /** \brief  Number of rows */
  size_type row_count() const ;

  /** \brief  Number of entries / end of row map */
  size_type entry_count() const ;

  /** \brief  Begining of entries for a given row */
  template< typename iType >
  size_type row_entry_begin( const iType & row ) const ;

  /** \brief  End of entries for a given row */
  template< typename iType >
  size_type row_entry_end(   const iType & row ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Has column mapping : entry -> column */
  bool has_column() const ;

  /** \brief  Column mapping : entry -> column */
  template< typename iType >
  size_type column( const iType & entry ) const ;

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

template< class CrsMapType , typename InputType >
inline
typename Impl::Factory< CrsMapType , InputType >::output_type
create_crsmap( const std::string & label , const InputType & input )
{
  return Impl::Factory< CrsMapType , InputType >::create( label , input );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Kokkos

#include <impl/Kokkos_CrsMap_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CRSMAP_HPP */

