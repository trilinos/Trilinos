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

#ifndef KOKKOS_DEVICECUDA_MULTIVECTORVIEW_HPP
#define KOKKOS_DEVICECUDA_MULTIVECTORVIEW_HPP

#include <Kokkos_MultiVectorView.hpp>
#include <Kokkos_DeviceCuda.hpp>

#include <DeviceCuda/Kokkos_DeviceCuda_ParallelDriver.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_DeepCopy.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include <impl/Kokkos_MultiVectorView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template < typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceCuda , DeviceCuda ,
                           true  /* Same memory space */ ,
                           false /* Not contiguous */ > {
public:
  typedef MultiVectorView< ValueType , DeviceCuda > multivector_type ;

  static
  void run( const multivector_type & dst , const multivector_type & src )
  {
    const DeviceCuda::size_type n = dst.m_dim[ multivector_type::RankStride ] *
                                    dst.m_dim[ multivector_type::RankCount ] ; 
    parallel_for( n , DeepCopyContiguous< ValueType , DeviceCuda >
                        ( dst.m_memory.ptr_on_device() ,
                          src.m_memory.ptr_on_device() ) );
  }
};

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceCuda , DeviceHost ,
                           false /* different memory space */ ,
                           false /* one is not contiguous */  >
{
public:
  typedef MultiVectorView< ValueType , DeviceCuda > dst_type ;
  typedef MultiVectorView< ValueType , DeviceHost > src_type ;

  inline
  static void run( const dst_type & dst , const src_type & src )
  {
    for ( DeviceCuda::size_type i = 0 ; i < dst.count() ; ++i ) {
      Impl::copy_to_cuda_from_host(
        dst.m_ptr_on_device + i * dst.m_dim[ dst_type::RankStride ] ,
        src.m_ptr_on_device + i * src.m_dim[ src_type::RankLength ] ,
        sizeof(ValueType), dst.length() );
    }
  }
};

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceHost , DeviceCuda ,
                           false /* different memory space */ ,
                           false /* one is not contiguous */  >
{
public:
  typedef MultiVectorView< ValueType , DeviceHost > dst_type ;
  typedef MultiVectorView< ValueType , DeviceCuda > src_type ;

  inline
  static void run( const dst_type & dst , const src_type & src )
  {
    for ( DeviceCuda::size_type i = 0 ; i < src.count() ; ++i ) {
      Impl::copy_to_host_from_cuda(
        dst.m_ptr_on_device + i * dst.m_dim[ dst_type::RankLength ] ,
        src.m_ptr_on_device + i * src.m_dim[ src_type::RankStride ] ,
        sizeof(ValueType), src.length() );
    }
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_DEVICECUDA_MULTIVECTORVIEW_HPP */


