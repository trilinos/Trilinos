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

#ifndef KOKKOS_DEVICE_CUDA_DEEP_COPY_HPP
#define KOKKOS_DEVICE_CUDA_DEEP_COPY_HPP

#include <Kokkos_ArrayForwardDeclarations.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>

namespace Kokkos {
namespace Impl {

void copy_to_cuda_from_host( void * dst , const void * src ,
                             size_t member_size , size_t member_count );

void copy_to_host_from_cuda( void * dst , const void * src ,
                             size_t member_size , size_t member_count );
}

/*------------------------------------------------------------------------*/

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceCuda >
{
public:
  static void run( const ValueView< ValueType , DeviceCuda > & dst ,
                   const ValueType & src )
  {
    ValueType * const d = dst.m_memory.ptr_on_device();
    Impl::copy_to_cuda_from_host( d , & src, sizeof(ValueType), 1 );
  }

  static void run( ValueType & dst ,
                   const ValueView< ValueType , DeviceCuda > & src )
  {
    ValueType * const s = src.m_memory.ptr_on_device();
    Impl::copy_to_host_from_cuda( & dst , s, sizeof(ValueType), 1 );
  }
};
/*------------------------------------------------------------------------*/

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceCuda , DeviceHost >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceCuda > & dst ,
                   const MultiVectorView< ValueType , DeviceHost > & src )
  {
    Impl::copy_to_cuda_from_host( dst.m_ptr_on_device ,
                                  src.m_ptr_on_device,
                                  sizeof(ValueType),
                                  dst.size() );
  }
};

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceHost , DeviceCuda >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceHost > & dst ,
                   const MultiVectorView< ValueType , DeviceCuda > & src )
  {
    Impl::copy_to_host_from_cuda( dst.m_ptr_on_device ,
                                  src.m_ptr_on_device,
                                  sizeof(ValueType),
                                  dst.size() );
  }
};

/*------------------------------------------------------------------------*/
/** \brief  Copy Host to Cuda specialization */
template< typename ValueType , class MapOpt >
class MDArrayDeepCopy< ValueType ,
                       DeviceCuda , MapOpt , true ,
                       DeviceHost , MapOpt , true >
{
public:
  typedef MDArrayView< ValueType , DeviceCuda , MapOpt > dst_type ;
  typedef MDArrayView< ValueType , DeviceHost , MapOpt > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    Impl::copy_to_cuda_from_host( dst.m_memory.ptr_on_device() ,
                                  src.m_memory.ptr_on_device() ,
                                  sizeof(ValueType) , dst.size() );
  }
};


/** \brief  Copy Cuda to Host specialization */
template< typename ValueType , class MapOpt >
class MDArrayDeepCopy< ValueType ,
                       DeviceHost , MapOpt , true ,
                       DeviceCuda , MapOpt , true >
{
public:
  typedef MDArrayView< ValueType , DeviceHost , MapOpt > dst_type ;
  typedef MDArrayView< ValueType , DeviceCuda , MapOpt > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    Impl::copy_to_host_from_cuda( dst.m_memory.ptr_on_device() ,
                                  src.m_memory.ptr_on_device() ,
                                  sizeof(ValueType) , dst.size() );
  }
};

/*------------------------------------------------------------------------*/

} // namespace Kokkos


#endif /* #ifndef KOKKOS_DEVICE_CUDA_DEEP_COPY_HPP */


