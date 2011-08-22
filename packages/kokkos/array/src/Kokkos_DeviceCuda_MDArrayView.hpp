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

#ifndef KOKKOS_DEVICECUDA_MDARRAYVIEW_HPP
#define KOKKOS_DEVICECUDA_MDARRAYVIEW_HPP

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_DeviceCuda.hpp>

#include <impl/Kokkos_ArrayBounds.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_ParallelDriver.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_DeepCopy.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapRight_macros.hpp>
#include <impl/Kokkos_MDArrayView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Deep copy device to device */
template< typename ValueType >
class MDArrayDeepCopy< ValueType , DeviceCuda , DeviceCuda ,
                       true  /* same memory space */ ,
                       true  /* same mdarray map */ ,
                       false /* contiguous memory */ >
{
public:
  typedef DeviceCuda            device_type ;
  typedef DeviceCuda::size_type size_type ;

  typedef MDArrayView< ValueType , DeviceCuda > array_type ;

        ValueType * const dst ;
  const ValueType * const src ;

  MDArrayDeepCopy( ValueType * arg_dst , const ValueType * arg_src )
    : dst( arg_dst ), src( arg_src ) {}

  inline
  __device__
  void operator()( size_type iwork ) const
  { dst[iwork] = src[iwork] ; }

  static void run( const array_type & dst , const array_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    const size_type n = dst.m_map.allocation_size();

    parallel_for( n , MDArrayDeepCopy( dst.m_memory.ptr_on_device() ,
                                       src.m_memory.ptr_on_device() ) );
  }
};

/** \brief  Copy Host to Cuda specialization */
template< typename ValueType >
class MDArrayDeepCopy< ValueType , DeviceCuda ,
                       Serial< HostMemory , MDArrayIndexMapLeft > ,
                       false /* same memory space */ ,
                       true  /* same mdarray map */ ,
                       false /* contiguous memory */ >
{
private:
  enum { OK = Impl::StaticAssert<
              Impl::SameType< MDArrayIndexMapLeft ,
                              DeviceCuda::mdarray_map >::value >::value };
public:
  typedef DeviceCuda::size_type                 size_type ;
  typedef MDArrayView< ValueType , DeviceCuda > dst_type ;
  typedef typename dst_type::HostView           src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    // The Cuda MDArray has MapLeft and the first dimension is padded.

    size_type count = 1 ;
    for ( size_type r = 1 ; r < dst.rank() ; ++r ) {
      count *= dst.m_map.dimension(r);
    }

    const size_type length = dst.m_map.dimension(0);
    const size_type stride = dst.m_map.allocation_size() / count ;

    for ( size_type i = 0 ; i < count ; ++i ) {
      Impl::copy_to_cuda_from_host(
        dst.m_memory.ptr_on_device() + stride * i ,
        src.m_memory.ptr_on_device() + length * i ,
        sizeof(ValueType) , length );
    }
  }
};


/** \brief  Copy Cuda to Host specialization */
template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       Serial< HostMemory , MDArrayIndexMapLeft > ,
                       DeviceCuda ,
                       false /* same memory space */ ,
                       true  /* same mdarray map */ ,
                       false /* contiguous memory */ >
{
private:
  enum { OK = Impl::StaticAssert<
              Impl::SameType< MDArrayIndexMapLeft ,
                              DeviceCuda::mdarray_map >::value >::value };
public:
  typedef DeviceCuda::size_type                 size_type ;
  typedef MDArrayView< ValueType , DeviceCuda > src_type ;
  typedef typename src_type::HostView           dst_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    // The Cuda MDArray has MapLeft and the first dimension is padded.

    size_type count = 1 ;
    for ( size_type r = 1 ; r < src.rank() ; ++r ) {
      count *= src.m_map.dimension(r);
    }

    const size_type length = src.m_map.dimension(0);
    const size_type stride = src.m_map.allocation_size() / count ;

    for ( size_type i = 0 ; i < count ; ++i ) {
      Impl::copy_to_host_from_cuda(
        dst.m_memory.ptr_on_device() + length * i ,
        src.m_memory.ptr_on_device() + stride * i ,
        sizeof(ValueType) , length );
    }
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_DEVICECUDA_MDARRAYVIEW_HPP */


