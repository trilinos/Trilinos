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

#ifndef KOKKOS_DEVICETPI_MDARRAYVIEW_HPP
#define KOKKOS_DEVICETPI_MDARRAYVIEW_HPP

#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_DeviceHost_MDArrayView.hpp>

#include <Kokkos_DeviceTPI_macros.hpp>
#include <impl/Kokkos_MDArrayView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Kokkos {
namespace Impl {

/*------------------------------------------------------------------------*/
/** \brief  Copy Host to TPI specialization with same map and contiguous */

template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       DeviceTPI  ,
                       Serial< HostMemory , DeviceTPI::mdarray_map > ,
                       true /* Same memory space */ ,
                       true /* Same map */ ,
                       true /* Contiguous */ >
{
private:
  typedef Serial< HostMemory , DeviceTPI::mdarray_map > device_host ;
public:
  typedef MDArrayView< ValueType , DeviceTPI   > dst_type ;
  typedef MDArrayView< ValueType , device_host > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef DeepCopyContiguous<ValueType,DeviceTPI> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};


/** \brief  Copy TPI to Host specialization with same map and contiguou */
template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       Serial< HostMemory , DeviceTPI::mdarray_map > ,
                       DeviceTPI ,
                       true /* Same memory space */ ,
                       true /* Same map */ ,
                       true /* Contiguous */ >
{
private:
  typedef Serial< HostMemory , DeviceTPI::mdarray_map > device_host ;
public:
  typedef MDArrayView< ValueType , device_host > dst_type ;
  typedef MDArrayView< ValueType , DeviceTPI >   src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef DeepCopyContiguous<ValueType,DeviceTPI> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

/*------------------------------------------------------------------------*/

} // namespace Impl
} // namespace Kokkos


#endif /* #ifndef KOKKOS_DEVICETPI_MDARRAYVIEW_HPP */

