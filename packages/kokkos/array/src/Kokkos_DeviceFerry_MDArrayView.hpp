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

#ifndef KOKKOS_DEVICEFERRY_MDARRAYVIEW_HPP
#define KOKKOS_DEVICEFERRY_MDARRAYVIEW_HPP

#include <Kokkos_MDArrayView.hpp>

#include <Kokkos_DeviceFerry_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapRight_macros.hpp>
#include <impl/Kokkos_MDArrayView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Kokkos {
namespace Impl {

/** \brief  Copy Host to Ferry specialization */
template< typename ValueType , class MapOpt >
class MDArrayDeepCopy< ValueType ,
                       DeviceFerry , MapOpt , true ,
                       DeviceHost , MapOpt , true >
{
public:
  typedef MDArrayView< ValueType , DeviceFerry , MapOpt > dst_type ;
  typedef MDArrayView< ValueType , DeviceHost , MapOpt > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );
    
	const ValueType * s = src.m_memory.ptr_on_device();
    long int d = (long int)dst.m_memory.ptr_on_device();
    int size = src.size();
    #pragma offload target(mic) in(s:length(size)) in(size) in(d)
    {
    	ValueType * local_d = (ValueType*)d;
    	for(int i = 0; i < size ; i++) {
    		local_d[i] = s[i];	
    	}
    }
/*    Impl::copy_to_ferry_from_host( dst.m_memory.ptr_on_device() ,
                                  src.m_memory.ptr_on_device() ,
                                  sizeof(ValueType) , dst.size() ); */
  }
};


/** \brief  Copy Ferry to Host specialization */
template< typename ValueType , class MapOpt >
class MDArrayDeepCopy< ValueType ,
                       DeviceHost , MapOpt , true ,
                       DeviceFerry , MapOpt , true >
{
public:
  typedef MDArrayView< ValueType , DeviceHost , MapOpt > dst_type ;
  typedef MDArrayView< ValueType , DeviceFerry , MapOpt > src_type ;

  static void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );
	long int s = (long int)src.m_memory.ptr_on_device();
    ValueType * d = dst.m_memory.ptr_on_device();
    int size = src.size();
    #pragma offload target(mic) in(s) in(size) out(d:length(size))
    {
    	ValueType * local_s = (ValueType*)s;
    	for(int i = 0; i < size ; i++) {
    		d[i] = local_s[i];	
    	}
    }
/*    Impl::copy_to_host_from_ferry( dst.m_memory.ptr_on_device() ,
                                  src.m_memory.ptr_on_device() ,
                                  sizeof(ValueType) , dst.size() ); */
  }
};

/*------------------------------------------------------------------------*/

} // namespace Impl
} // namespace Kokkos


#endif /* #ifndef KOKKOS_DEVICEFERRY_MDARRAYVIEW_HPP */


