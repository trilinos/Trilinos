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

#ifndef KOKKOS_HOSTMAP_HPP
#define KOKKOS_HOSTMAP_HPP

#include <cstddef>
#include <list>
#include <Kokkos_BaseMappedArray.hpp>

namespace Kokkos {

class HostDevice ;
class HostMap ;

/*--------------------------------------------------------------------------*/

class HostDevice : public BaseDeviceMemory {
public:
  typedef HostMap map_type ;

  HostDevice();

  void * allocate( size_type value_size ,
                   size_type chunk_count ,
                   size_type work_count );

  void deallocate( void * );

  ~HostDevice();

private:
  HostDevice( const HostDevice & );
  HostDevice operator = ( const HostDevice & );
};

/*--------------------------------------------------------------------------*/

class HostMap {
public:
  typedef HostDevice                 device_type ;
  typedef BaseMappedArray::size_type size_type ;

  /*------------------------------------------------------------------------*/

  template< typename Scalar ,
            typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 , const iType1 & i1 ,
                      const iType2 & i2 , const iType3 & i3 ,
                      const iType4 & i4 , const iType5 & i5 ,
                      const iType6 & i6 , const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds(i0,i1,i2,i3,i4,i5,i6,iP) );

      Scalar          * const ptr = array.pointer_on_device<Scalar>();
      const size_type * const dim = array.dimension();

      return ptr[ ( i0 + dim[0] * ( i1 + dim[1] *
                  ( i2 + dim[2] * ( i3 + dim[3] *
                  ( i4 + dim[4] * ( i5 + dim[5] *
                  ( i6 + dim[6] * ( iP )))))))) ];
    }

  template< typename Scalar ,
            typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 , const iType1 & i1 ,
                      const iType2 & i2 , const iType3 & i3 ,
                      const iType4 & i4 , const iType5 & i5 ,
                      const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds(i0,i1,i2,i3,i4,i5,iP) );

      Scalar          * const ptr = array.pointer_on_device<Scalar>();
      const size_type * const dim = array.dimension();

      return ptr[ ( i0 + dim[0] * ( i1 + dim[1] *
                  ( i2 + dim[2] * ( i3 + dim[3] *
                  ( i4 + dim[4] * ( i5 + dim[5] *
                  ( iP ))))))) ];
    }

  template< typename Scalar ,
            typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 , const iType1 & i1 ,
                      const iType2 & i2 , const iType3 & i3 ,
                      const iType4 & i4 , const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds(i0,i1,i2,i3,i4,iP) );

      Scalar          * const ptr = array.pointer_on_device<Scalar>();
      const size_type * const dim = array.dimension();

      return ptr[ ( i0 + dim[0] * ( i1 + dim[1] *
                  ( i2 + dim[2] * ( i3 + dim[3] *
                  ( i4 + dim[4] * ( iP )))))) ];
    }

  template< typename Scalar ,
            typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 , const iType1 & i1 ,
                      const iType2 & i2 , const iType3 & i3 ,
                      const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds(i0,i1,i2,i3,iP) );

      Scalar          * const ptr = array.pointer_on_device<Scalar>();
      const size_type * const dim = array.dimension();

      return ptr[ ( i0 + dim[0] * ( i1 + dim[1] *
                  ( i2 + dim[2] * ( i3 + dim[3] *
                  ( iP ))))) ];
    }

  template< typename Scalar ,
            typename iType0 , typename iType1 ,
            typename iType2 , typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 , const iType1 & i1 ,
                      const iType2 & i2 , const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds(i0,i1,i2,iP) );

      Scalar          * const ptr = array.pointer_on_device<Scalar>();
      const size_type * const dim = array.dimension();

      return ptr[ ( i0 + dim[0] * ( i1 + dim[1] *
                  ( i2 + dim[2] * ( iP )))) ];
    }

  template< typename Scalar ,
            typename iType0 , typename iType1 ,
            typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 , const iType1 & i1 ,
                      const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds(i0,i1,iP) );

      Scalar          * const ptr = array.pointer_on_device<Scalar>();
      const size_type * const dim = array.dimension();

      return ptr[ ( i0 + dim[0] * ( i1 + dim[1] * ( iP ))) ];
    }

  template< typename Scalar ,
            typename iType0 , typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iType0 & i0 , const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK(
        array.require_multi_index_in_bounds(i0,iP) );

      Scalar          * const ptr = array.pointer_on_device<Scalar>();
      const size_type * const dim = array.dimension();

      return ptr[ ( i0 + dim[0] * ( iP )) ];
    }

  template< typename Scalar , typename iTypeP >
  inline static
  Scalar & value_map( const BaseMappedArray & array ,
                      const iTypeP & iP )
    {
      KOKKOS_MDARRAY_CHECK( array.require_multi_index_in_bounds(iP) );
      return array.pointer_on_device<Scalar>()[ iP ];
    }

  /*------------------------------------------------------------------------*/
  /** \brief  Map arrays of parallel_work_length onto device.  */
  HostMap( size_type parallel_work_count );

  /** \brief  Destroy the host map and all of its allocated arrays. */
  ~HostMap();

  inline
  void allocate( BaseMappedArray & array ,
                 size_type sizeof_value ,
                 size_type rank ,
                 const size_type dimension[] )
    { m_arrays.allocate( array , sizeof_value , rank , dimension ); }


  inline
  void deallocate( BaseMappedArray & array )
    { m_arrays.deallocate( array ); }

private:

  HostDevice                m_device ; ///< Device for allocated array memory
  BaseMappedArrayRepository m_arrays ; ///< Repository depends on device

  HostMap();
  HostMap( const HostMap & );
  HostMap & operator = ( const HostMap & );
};

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

#endif /* KOKKOS_HOSTMAP_HPP */


