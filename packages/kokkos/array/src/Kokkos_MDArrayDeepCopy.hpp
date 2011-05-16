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

#ifndef KOKKOS_MDARRAYDEEPCOPY_HPP
#define KOKKOS_MDARRAYDEEPCOPY_HPP

#include <Kokkos_MDArrayView.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

namespace Kokkos {

template< typename ValueType ,
          class DeviceDst , class MapDst ,
          class DeviceSrc , class MapSrc >
class MDArrayDeepCopy {
private:
  enum { DifferentDevices = ! Impl::SameType<DeviceDst,DeviceSrc>::value };
  enum { DifferentMaps    = ! Impl::SameType<MapDst,MapSrc>::value };
  enum { OK = Impl::StaticAssert< DifferentDevices && DifferentMaps >::OK };
public:

  typedef MDArrayView<ValueType,DeviceDst,MapDst> dst_type ;
  typedef MDArrayView<ValueType,DeviceSrc,MapSrc> src_type ;

  // Both the devices and the maps are different.
  // Have to make a temporary copy to bring into alignment.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    // Copy to the destination device with the same map:
    MDArrayView<ValueType,DeviceDst,MapSrc> tmp_src ;
    
    MDArrayDeepCopy<ValueType,DeviceDst,MapSrc,
                              DeviceSrc,MapSrc>
      ::run( tmp_src , src );

    // Copy with remapping on the destination device:
    MDArrayDeepCopy<ValueType,DeviceDst,MapDst,
                              DeviceDst,MapSrc>
      ::run( dst , tmp_src );
  }
};

template< typename ValueType ,
          class DeviceDst , class MapDst ,
          class DeviceSrc , class MapSrc >
void deep_copy( const MDArrayView<ValueType,DeviceDst,MapDst> & dst ,
                const MDArrayView<ValueType,DeviceSrc,MapSrc> & src )
{
  MDArrayDeepCopy<ValueType,DeviceDst,MapDst,DeviceSrc,MapSrc>
    ::run( dst , src );
}

} // namespace Kokkos

//----------------------------------------------------------------------------

#if defined( KOKKOS_DEVICE_HOST )
#include <DeviceHost/Kokkos_HostMDArrayDeepCopy.hpp>
#endif

//----------------------------------------------------------------------------

#endif /* KOKKOS_MDARRAYDEEPCOPY_HPP */


