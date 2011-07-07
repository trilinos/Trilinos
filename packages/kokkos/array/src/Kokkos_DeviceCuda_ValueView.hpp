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

#ifndef KOKKOS_DEVICE_CUDA_VALUEVIEW_HPP
#define KOKKOS_DEVICE_CUDA_VALUEVIEW_HPP

#include <Kokkos_ValueView.hpp>

#include <Kokkos_DeviceCuda_macros.hpp>
#include <impl/Kokkos_ValueView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

#include <DeviceCuda/Kokkos_DeviceCuda_DeepCopy.hpp>

namespace Kokkos {
namespace Impl {

/*------------------------------------------------------------------------*/

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceCuda , DeviceHost >
{
public:
  static void run( const ValueView< ValueType , DeviceCuda > & dst ,
                   const ValueType & src )
  {
    ValueType * const d = dst.m_memory.ptr_on_device();
    Impl::copy_to_cuda_from_host( d , & src, sizeof(ValueType), 1 );
  }

  static void run( const ValueView< ValueType , DeviceCuda > & dst ,
                   const ValueView< ValueType , DeviceHost > & src )
  {
    ValueType * const d = dst.m_memory.ptr_on_device();
    ValueType * const s = src.m_memory.ptr_on_device();
    Impl::copy_to_cuda_from_host( d , s, sizeof(ValueType), 1 );
  }
};

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceHost , DeviceCuda >
{
public:

  static void run( ValueType & dst ,
                   const ValueView< ValueType , DeviceCuda > & src )
  {
    ValueType * const s = src.m_memory.ptr_on_device();
    Impl::copy_to_host_from_cuda( & dst , s, sizeof(ValueType), 1 );
  }

  static void run( const ValueView< ValueType , DeviceHost > & dst ,
                   const ValueView< ValueType , DeviceCuda > & src )
  {
    ValueType * const d = dst.m_memory.ptr_on_device();
    ValueType * const s = src.m_memory.ptr_on_device();
    Impl::copy_to_host_from_cuda( d , s, sizeof(ValueType), 1 );
  }
};

/*------------------------------------------------------------------------*/

} // namespace Impl
} // namespace Kokkos


#endif /* #ifndef KOKKOS_DEVICECUDA_VALUEVIEW_HPP */

