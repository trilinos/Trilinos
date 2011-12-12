/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef KOKKOS_DEVICEPTHREAD_MDARRAYVIEW_HPP
#define KOKKOS_DEVICEPTHREAD_MDARRAYVIEW_HPP

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>

#include <Kokkos_DevicePthread_macros.hpp>
#include <impl/Kokkos_MDArrayView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Kokkos {
namespace Impl {

/*------------------------------------------------------------------------*/
/** \brief  Copy Host to Pthread specialization with same map and contiguous */

template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       DevicePthread  ,
                       Serial< HostMemory , DevicePthread::mdarray_map > ,
                       true /* Same memory space */ ,
                       true /* Same map */ ,
                       true /* Contiguous */ >
{
private:
  typedef Serial< HostMemory , DevicePthread::mdarray_map > device_host ;
public:
  typedef MDArrayView< ValueType , DevicePthread   > dst_type ;
  typedef MDArrayView< ValueType , device_host > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef DeepCopyContiguous<ValueType,DevicePthread> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};


/** \brief  Copy Pthread to Host specialization with same map and contiguou */
template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       Serial< HostMemory , DevicePthread::mdarray_map > ,
                       DevicePthread ,
                       true /* Same memory space */ ,
                       true /* Same map */ ,
                       true /* Contiguous */ >
{
private:
  typedef Serial< HostMemory , DevicePthread::mdarray_map > device_host ;
public:
  typedef MDArrayView< ValueType , device_host > dst_type ;
  typedef MDArrayView< ValueType , DevicePthread >   src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef DeepCopyContiguous<ValueType,DevicePthread> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

/*------------------------------------------------------------------------*/

} // namespace Impl
} // namespace Kokkos


#endif /* #ifndef KOKKOS_DEVICEPTHREAD_MDARRAYVIEW_HPP */

