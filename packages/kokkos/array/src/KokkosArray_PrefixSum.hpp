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

#ifndef KOKKOS_PREFIXSUM_HPP
#define KOKKOS_PREFIXSUM_HPP

#include <string>
#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_ArrayBounds.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------
/** \brief  Prefix sum of integer values.  */

template< typename IntType , class DeviceType >
class PrefixSum {
public:
  typedef DeviceType  device_type ;
  typedef IntType     size_type ;

  typedef PrefixSum< size_type , typename HostMapped< device_type >::type >
          HostMirror ;

  /*------------------------------------------------------------------*/
  /** \brief  Number of values  */
  size_type length() const ;

  /** \brief  Sum of values == operator[]( count() ) */
  size_type sum() const ;

  /** \brief  Begining of entries for a given row */
  template< typename iType >
  size_type operator[]( const iType & ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  PrefixSum();

  /** \brief  Construct a view of the array */
  PrefixSum( const PrefixSum & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  PrefixSum & operator = ( const PrefixSum & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~PrefixSum();

  /*------------------------------------------------------------------*/
  operator bool() const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const PrefixSum & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const PrefixSum & ) const ;
};

//----------------------------------------------------------------------------

template< class PrefixSumType , class InputType >
inline
PrefixSum< typename PrefixSumType::size_type ,
           typename PrefixSumType::device_type >
create_prefixsum( const std::string & label ,
                  const InputType & input )
{
  return Impl::Factory< PrefixSumType , InputType >::create( label , input );
}

template< class PrefixSumType , class InputType >
inline
PrefixSum< typename PrefixSumType::size_type ,
           typename PrefixSumType::device_type >
create_prefixsum( const InputType & input )
{
  return Impl::Factory< PrefixSumType , InputType >
             ::create( std::string() , input );
}

//----------------------------------------------------------------------------

template< typename IntType ,
          class DeviceDst ,
          class DeviceSrc >
inline
void deep_copy(       PrefixSum<IntType,DeviceDst> & dst ,
                const PrefixSum<IntType,DeviceSrc> & src )
{
  typedef PrefixSum<IntType,DeviceDst> dst_type ;
  typedef PrefixSum<IntType,DeviceSrc> src_type ;

  if ( dst.operator!=(src) ) {

    Impl::prefixsum_require_equal_dimension( dst.length() , src.length() );

    Impl::Factory< dst_type , src_type >::deep_copy( dst , src );
  }
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/KokkosArray_PrefixSum_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_PREFIXSUM_HPP */

