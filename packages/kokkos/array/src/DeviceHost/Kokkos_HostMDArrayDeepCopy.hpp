/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_HOSTMDARRAYDEEPCOPY_HPP
#define KOKKOS_HOSTMDARRAYDEEPCOPY_HPP

#include <Kokkos_MDArrayDeepCopy.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

#include <impl/Kokkos_DeviceHost_macros.hpp>
#include <impl/Kokkos_MDArrayDeepCopy_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>

namespace Kokkos {

template< typename ValueType , class MapType >
class MDArrayDeepCopy< ValueType , DeviceHost , MapType ,
                                   DeviceHost , MapType >
{
public:
  typedef MDArrayView< ValueType , DeviceHost , MapType > mdarray_type ;

  static
  void run( const mdarray_type & dst , const mdarray_type & src )
  {
    if ( mdarray_type::Contiguous ) {
      size_t dst_dims[8] ;
      size_t src_dims[8] ;
      dst.dimensions( dst_dims );
      src.dimensions( src_dims );

      Impl::mdarray_require_equal_dimension(
        dst.rank() , dst_dims , src.rank() , src_dims );

      const ValueType * src_ptr = src.ptr_on_device();
            ValueType * dst_ptr = dst.ptr_on_device();
            ValueType * const dst_end = dst_ptr + dst.size();
      while ( dst_ptr < dst_end ) { *dst_ptr++ = *src_ptr++ ; }
    }
    else {

    }
  }
};

template< typename ValueType , class MapDst , class MapSrc >
class MDArrayDeepCopy< ValueType , DeviceHost , MapDst ,
                                   DeviceHost , MapSrc >
{
public:
  typedef DeviceHost::size_type size_type ;

  typedef MDArrayView< ValueType , DeviceHost , MapSrc > src_type ;
  typedef MDArrayView< ValueType , DeviceHost , MapDst > dst_type ;

  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 8 > deep8 ;
  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 7 > deep7 ;
  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 6 > deep6 ;
  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 5 > deep5 ;
  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 4 > deep4 ;
  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 3 > deep3 ;
  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 2 > deep2 ;
  Impl::MDArrayDeepCopyMember<ValueType,DeviceHost, MapDst, MapSrc, 1 > deep1 ;

  static
  void run( const dst_type & dst , const src_type & src )
  {
    size_t dst_dims[8] ;
    size_t src_dims[8] ;
    dst.dimensions( dst_dims );
    src.dimensions( src_dims );
    Impl::mdarray_require_equal_dimension(
      dst.rank() , dst_dims , src.rank() , src_dims );


    const size_type n = dst.size();

    switch ( dst.rank() ) {
    case 8 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep8:copy( i , dst , src ); }
      break ;
    case 7 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep7:copy( i , dst , src ); }
      break ;
    case 6 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep6:copy( i , dst , src ); }
      break ;
    case 5 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep5:copy( i , dst , src ); }
      break ;
    case 4 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep4:copy( i , dst , src ); }
      break ;
    case 3 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep3:copy( i , dst , src ); }
      break ;
    case 2 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep2:copy( i , dst , src ); }
      break ;
    case 1 :
      for ( size_type i = 0 ; i < n ; ++i ) { deep1:copy( i , dst , src ); }
      break ;
    }
  }
};

} // namespace Kokkos

#endif /* KOKKOS_HOSTMDARRAYDEEPCOPY_HPP */


