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

#include <Kokkos_DeviceHost_ParallelFor.hpp>

namespace Kokkos {
namespace Impl {
  
template< class DeviceType , class MDArrayDst ,
                             class MDArraySrc , unsigned Rank >
class MDArrayDeepCopyFunctor ;

//----------------------------------------------------------------------------

/** \brief  Fully conformal and contiguous arrays */
template< typename ValueType , class DeviceType >
class MDArrayDeepCopy< ValueType , DeviceType , DeviceType ,
                       true , true , true >
{
public:

  typedef MDArrayView<ValueType,DeviceType> array_type ;

  static void run( const array_type & dst , const array_type & src )
  {
    typedef Impl::DeepCopyContiguous<ValueType,DeviceType> functor_type ;

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

//----------------------------------------------------------------------------
/** \brief  Deep copy from Device to Host,
 *          with different maps,
 *          and no assumption of contiguity.
 *          Force remap to occur on the host.
 */
template< typename ValueType , class DeviceSrc , bool Contig >
class MDArrayDeepCopy< ValueType , DeviceHost , DeviceSrc ,
                       false  /* Different Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
private:
  typedef typename DeviceHost::mdarray_map  map_dst ;
  typedef typename DeviceSrc ::mdarray_map  map_src ;
public:

  typedef Serial< HostMemory , map_src >  DeviceTmp ;

  typedef MDArrayView<ValueType,DeviceHost> dst_type ;
  typedef MDArrayView<ValueType,DeviceTmp>  tmp_type ;
  typedef MDArrayView<ValueType,DeviceSrc>  src_type ;

  // Both the devices and the maps are different.
  // Copy to a temporary on the destination with the source map
  // and then remap the temporary to the final array.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    size_t dims[ MDArrayMaxRank ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    dst.dimensions( dims );

    tmp_type tmp = create_labeled_mdarray<tmp_type>(
                                "temporary" ,
                                dims[0] , dims[1] , dims[2] , dims[3] ,
                                dims[4] , dims[5] , dims[6] , dims[7] );

    MDArrayDeepCopy< ValueType , DeviceTmp ,  DeviceSrc >::run( tmp , src );
    MDArrayDeepCopy< ValueType , DeviceHost , DeviceTmp >::run( dst , tmp );
  }
};

template< typename ValueType , class MapDst , class DeviceSrc , bool Contig >
class MDArrayDeepCopy< ValueType , Serial< HostMemory , MapDst > , DeviceSrc ,
                       false  /* Different Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
private:

  typedef Serial< HostMemory , MapDst >     device_host ;
  typedef typename DeviceSrc ::mdarray_map  map_src ;
public:

  typedef Serial< HostMemory , map_src >  DeviceTmp ;

  typedef MDArrayView<ValueType,device_host> dst_type ;
  typedef MDArrayView<ValueType,DeviceTmp>   tmp_type ;
  typedef MDArrayView<ValueType,DeviceSrc>   src_type ;

  // Both the devices and the maps are different.
  // Copy to a temporary on the destination with the source map
  // and then remap the temporary to the final array.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    size_t dims[ MDArrayMaxRank ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    dst.dimensions( dims );

    tmp_type tmp = create_labeled_mdarray<tmp_type>(
                                "temporary" ,
                                dims[0] , dims[1] , dims[2] , dims[3] ,
                                dims[4] , dims[5] , dims[6] , dims[7] );

    MDArrayDeepCopy< ValueType , DeviceTmp ,  DeviceSrc >::run( tmp , src );
    MDArrayDeepCopy< ValueType , device_host, DeviceTmp >::run( dst , tmp );
  }
};

//----------------------------------------------------------------------------
/** \brief  Deep copy from Host to Device,
 *          with different maps,
 *          and no assumption of contiguity.
 *          Force remap to occur on the host.
 */
template< typename ValueType , class DeviceDst , bool Contig >
class MDArrayDeepCopy< ValueType , DeviceDst , DeviceHost ,
                       false  /* Different Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
private:
  typedef typename DeviceDst ::mdarray_map  map_dst ;
  typedef typename DeviceHost::mdarray_map  map_src ;
public:

  typedef Serial< HostMemory , map_dst >  DeviceTmp ;

  typedef MDArrayView<ValueType,DeviceDst>  dst_type ;
  typedef MDArrayView<ValueType,DeviceTmp>  tmp_type ;
  typedef MDArrayView<ValueType,DeviceHost> src_type ;

  // Both the devices and the maps are different.
  // Copy to a temporary on the destination with the source map
  // and then remap the temporary to the final array.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    size_t dims[ MDArrayMaxRank ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    dst.dimensions( dims );

    tmp_type tmp = create_labeled_mdarray<tmp_type>(
                                "temporary" ,
                                dims[0] , dims[1] , dims[2] , dims[3] ,
                                dims[4] , dims[5] , dims[6] , dims[7] );

    MDArrayDeepCopy< ValueType , DeviceTmp , DeviceHost >::run( tmp , src );
    MDArrayDeepCopy< ValueType , DeviceDst , DeviceTmp > ::run( dst , tmp );
  }
};

template< typename ValueType , class DeviceDst , class MapSrc , bool Contig >
class MDArrayDeepCopy< ValueType , DeviceDst , Serial< HostMemory , MapSrc > ,
                       false  /* Different Memory Space */ ,
                       false  /* Different MDArray Maps */ ,
                       Contig /* Don't care */ >
{
private:
  typedef typename DeviceDst ::mdarray_map  map_dst ;
  typedef Serial< HostMemory , MapSrc >     device_host ;
public:

  typedef Serial< HostMemory , map_dst >  DeviceTmp ;

  typedef MDArrayView<ValueType,DeviceDst>   dst_type ;
  typedef MDArrayView<ValueType,DeviceTmp>   tmp_type ;
  typedef MDArrayView<ValueType,device_host> src_type ;

  // Both the devices and the maps are different.
  // Copy to a temporary on the destination with the source map
  // and then remap the temporary to the final array.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    size_t dims[ MDArrayMaxRank ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    dst.dimensions( dims );

    tmp_type tmp = create_labeled_mdarray<tmp_type>(
                                "temporary" ,
                                dims[0] , dims[1] , dims[2] , dims[3] ,
                                dims[4] , dims[5] , dims[6] , dims[7] );

    MDArrayDeepCopy< ValueType , DeviceTmp , device_host >::run( tmp , src );
    MDArrayDeepCopy< ValueType , DeviceDst , DeviceTmp  > ::run( dst , tmp );
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_MDARRAYDEEPCOPY_HPP */


