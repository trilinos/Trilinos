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
    if ( MapType::contiguous ) {
      Impl::require_equal_dimension( dst , src );
      const ValueType * src = src.ptr_on_device();
            ValueType * dst = dst.ptr_on_device();
            ValueType * const dst_end = dst + dst.size();
      while ( dst < dst_end ) { *dst++ = *src++ ; }
    }
  }
};

template< typename ValueType , class MapTypeDst , class MapTypeSrc >
class MDArrayDeepCopy< ValueType , DeviceHost , MapTypeDst ,
                                   DeviceHost , MapTypeSrc >
{
public:
  typedef MDArrayView< ValueType , DeviceHost , MapTypeSrc > src_type ;
  typedef MDArrayView< ValueType , DeviceHost , MapTypeDst > dst_type ;

  static
  void run( const mdarray_type & dst , const mdarray_type & src )
  {
    if ( MapType::composition ) {
      Impl::require_equal_dimension( dst , src );

      const DeviceHost::size_type n = dst.size();
      DeviceHost::size_type indices[ dst_type::MAX_RANK ];

      ValueType * const dst = dst.ptr_on_device();

      switch( dst.rank() ) {
      case 8 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          dst.inverse_map( i , indices );

          dst( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] , indices[5] , indices[6] , indices[7] ) =
          src( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] , indices[5] , indices[6] , indices[7] );
        }
        break ;

      case 7 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          dst.inverse_map( i , indices );

          dst( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] , indices[5] , indices[6] ) =
          src( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] , indices[5] , indices[6] );
        }
        break ;

      case 6 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          dst.inverse_map( i , indices );

          dst( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] , indices[5] ) =
          src( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] , indices[5] );
        }
        break ;

      case 5 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          dst.inverse_map( i , indices );

          dst( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] ) =
          src( indices[0] , indices[1] , indices[2] , indices[3] ,
               indices[4] );
        }
        break ;

      case 4 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          dst.inverse_map( i , indices );

          dst( indices[0] , indices[1] , indices[2] , indices[3] ) =
          src( indices[0] , indices[1] , indices[2] , indices[3] );
        }
        break ;

      case 3 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          dst.inverse_map( i , indices );

          dst( indices[0] , indices[1] , indices[2] ) =
          src( indices[0] , indices[1] , indices[2] );
        }
        break ;

      case 2 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          dst.inverse_map( i , indices );

          dst( indices[0] , indices[1] ) =
          src( indices[0] , indices[1] );
        }
        break ;

      case 1 :
        for ( DeviceHost::size_type i = 0 ; i < n ; ++i ) {
          src.inverse_map( i , indices );
          dst( indices[0] ) = src( indices[0] );
        }
        break ;
    }
  }
};

} // namespace Kokkos

#endif /* KOKKOS_HOSTMDARRAYDEEPCOPY_HPP */


