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

#ifndef KOKKOS_ARRAY_HPP
#define KOKKOS_ARRAY_HPP

#include <string>
#include <KokkosArray_View.hpp>
#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_ArrayBounds.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Compile-time dimensioned array.  */
template< class ArrayType , class DeviceType >
class Array {
public:
  typedef DeviceType  device_type ;
  typedef ArrayType   array_type ;
  typedef typename device_type::size_type size_type ;

  typedef typename Impl::remove_all_extents<array_type>::type  value_type ;

  typedef Array< array_type , typename HostMapped< device_type >::type >
          HostMirror ;

  static const unsigned Rank = Impl::rank<array_type>::value + 1 ;

  template< unsigned I >
  struct Dimension {
    static const unsigned value = Impl::extent<array_type,I-1>::value ;
  };

  inline
  size_type rank() const ;

  template< typename iType >
  inline
  size_type dimension( const iType & rank ) const ;

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
  Array();

  /** \brief  Construct a view of the array */
  Array( const Array & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  Array & operator = ( const Array & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~Array();

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const Array & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const Array & ) const ;
};

//----------------------------------------------------------------------------

template< class ArrayType >
inline
Array< typename ArrayType::array_type ,
       typename ArrayType::device_type >
create_array( const std::string & label , const size_t nP )
{
  return Impl::Factory< ArrayType , void >::create( label , nP );
}

template< class ArrayType >
inline
Array< typename ArrayType::array_type ,
       typename ArrayType::device_type >
create_array( const size_t nP )
{
  return Impl::Factory< ArrayType , void >::create( std::string() , nP );
}

//----------------------------------------------------------------------------

template< class ArrayType ,
          class DeviceDst ,
          class DeviceSrc >
inline
void deep_copy( const Array<ArrayType,DeviceDst> & dst ,
                const Array<ArrayType,DeviceSrc> & src )
{
  typedef Array<ArrayType,DeviceDst> dst_type ;
  typedef Array<ArrayType,DeviceSrc> src_type ;

  if ( dst.operator!=(src) ) {

    Impl::array_require_equal_dimension( dst.dimension(0) , src.dimension(0) );

    Impl::Factory< dst_type , src_type >::deep_copy( dst , src );
  }
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/KokkosArray_Array_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_ARRAY_HPP */

