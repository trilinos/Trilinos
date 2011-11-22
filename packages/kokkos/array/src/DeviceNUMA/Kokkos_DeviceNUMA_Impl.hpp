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

//----------------------------------------------------------------------------
// Interfaces always included

#include <DeviceNUMA/Kokkos_DeviceNUMA_Parallel.hpp>
#include <DeviceNUMA/Kokkos_DeviceNUMA_For.hpp>
#include <DeviceNUMA/Kokkos_DeviceNUMA_Reduce.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Partial specializations for the device

#include <Kokkos_DeviceNUMA_macros.hpp>

#if ! defined( KOKKOS_DEVICENUMA_BASICFUNCTORS )
#define KOKKOS_DEVICENUMA_BASICFUNCTORS
#include <impl/Kokkos_BasicFunctors_macros.hpp>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_VALUEVIEW_HPP ) && \
  ! defined( KOKKOS_DEVICENUMA_VALUEVIEW )

#define KOKKOS_DEVICENUMA_VALUEVIEW

#include <impl/Kokkos_ValueView_macros.hpp>

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceNUMA , DeviceHost > {
public:

  static void run( const ValueView< ValueType , DeviceNUMA > & dst ,
                   const ValueType & src )
  { *dst = src ; }

  static void run( const ValueView< ValueType , DeviceNUMA >  & dst ,
                   const ValueView< ValueType , DeviceHost > & src )
  { *dst = *src ; }
};

template< typename ValueType >
class ValueDeepCopy< ValueType , DeviceHost , DeviceNUMA > {
public:

  static void run( ValueType & dst ,
                   const ValueView< ValueType , DeviceNUMA >  & src )
  { dst = *src ; }

  static void run( const ValueView< ValueType , DeviceHost > & dst ,
                   const ValueView< ValueType , DeviceNUMA >  & src )
  { *dst = *src ; }
};

} // namespace Impl
} // namespace Kokkos

#endif /* #if defined( KOKKOS_VALUEVIEW_HPP ) && ! defined( KOKKOS_DEVICENUMA_VALUEVIEW ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_MULTIVECTORVIEW_HPP ) && \
  ! defined( KOKKOS_DEVICENUMA_MULTIVECTORVIEW )

#define KOKKOS_DEVICENUMA_MULTIVECTORVIEW

#include <impl/Kokkos_MultiVectorView_macros.hpp>

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceNUMA , DeviceHost ,
                           true /* same memory space */ ,
                           true /* both are contiguous */ >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceNUMA >  & dst ,
                   const MultiVectorView< ValueType , DeviceHost > & src )
  {
    parallel_for( dst.length() * dst.count() ,
                  DeepCopyContiguous< ValueType , DeviceNUMA >
                    ( dst.m_memory.ptr_on_device() ,
                      src.m_memory.ptr_on_device() ) );
  }
};

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceHost , DeviceNUMA ,
                           true /* same memory space */ ,
                           true /* both are contiguous */ >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceHost > & dst ,
                   const MultiVectorView< ValueType , DeviceNUMA >  & src )
  {
    parallel_for( dst.length() * dst.count() ,
                  DeepCopyContiguous< ValueType , DeviceNUMA >
                    ( dst.m_memory.ptr_on_device() ,
                      src.m_memory.ptr_on_device() ) );
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* #if defined( KOKKOS_MULTIVECTORVIEW_HPP ) && ! defined( KOKKOS_DEVICENUMA_MULTIVECTORVIEW ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_MDARRAYVIEW_HPP ) && \
  ! defined( KOKKOS_DEVICENUMA_MDARRAYVIEW )

#define KOKKOS_DEVICENUMA_MDARRAYVIEW

#include <impl/Kokkos_MDArrayView_macros.hpp>

namespace Kokkos {
namespace Impl {

/*------------------------------------------------------------------------*/
/** \brief  Copy Host to NUMA specialization with same map and contiguous */

template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       DeviceNUMA  ,
                       Serial< HostMemory , DeviceNUMA::mdarray_map > ,
                       true /* Same memory space */ ,
                       true /* Same map */ ,
                       true /* Contiguous */ >
{
private:
  typedef Serial< HostMemory , DeviceNUMA::mdarray_map > device_host ;
public:
  typedef MDArrayView< ValueType , DeviceNUMA   > dst_type ;
  typedef MDArrayView< ValueType , device_host > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef DeepCopyContiguous<ValueType,DeviceNUMA> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};


/** \brief  Copy NUMA to Host specialization with same map and contiguou */
template< typename ValueType >
class MDArrayDeepCopy< ValueType ,
                       Serial< HostMemory , DeviceNUMA::mdarray_map > ,
                       DeviceNUMA ,
                       true /* Same memory space */ ,
                       true /* Same map */ ,
                       true /* Contiguous */ >
{
private:
  typedef Serial< HostMemory , DeviceNUMA::mdarray_map > device_host ;
public:
  typedef MDArrayView< ValueType , device_host > dst_type ;
  typedef MDArrayView< ValueType , DeviceNUMA >   src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    typedef DeepCopyContiguous<ValueType,DeviceNUMA> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

/*------------------------------------------------------------------------*/

} // namespace Impl
} // namespace Kokkos

#endif /* #if defined( KOKKOS_MDARRAYVIEW_HPP ) && ! defined( KOKKOS_DEVICENUMA_MDARRAYVIEW ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <Kokkos_DeviceClear_macros.hpp>




