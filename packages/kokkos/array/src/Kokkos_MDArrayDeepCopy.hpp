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

#include <Kokkos_ArrayForwardDeclarations.hpp>
#include <Kokkos_MDArrayView.hpp>
#include <Kokkos_ParallelFor.hpp>

#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>


namespace Kokkos {

//----------------------------------------------------------------------------

template< typename ValueType ,
          class DeviceDst , class MapDst ,
          class DeviceSrc , class MapSrc >
inline
void deep_copy( const MDArrayView<ValueType,DeviceDst,MapDst> & dst ,
                const MDArrayView<ValueType,DeviceSrc,MapSrc> & src )
{
  typedef MDArrayView<ValueType,DeviceDst,MapDst>  dst_type ;
  typedef MDArrayView<ValueType,DeviceSrc,MapSrc>  src_type ;

  typedef MDArrayDeepCopy<ValueType,DeviceDst,MapDst,dst_type::Contiguous,
                                    DeviceSrc,MapSrc,src_type::Contiguous>
    deep_copy_operator ;

  deep_copy_operator::run( dst , src );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_DEVICE_HOST )
#include <impl/Kokkos_DeviceHost_macros.hpp>
#include <impl/Kokkos_MDArrayDeepCopy_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_TPI )
#include <DeviceTPI/Kokkos_DeviceTPI_DeepCopy.hpp>
#include <impl/Kokkos_DeviceTPI_macros.hpp>
#include <impl/Kokkos_MDArrayDeepCopy_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_CUDA )
#include <DeviceCuda/Kokkos_DeviceCuda_DeepCopy.hpp>
#include <impl/Kokkos_DeviceCuda_macros.hpp>
#include <impl/Kokkos_MDArrayDeepCopy_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>
#endif


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Deep copy between different devices, different maps,
 *          and no assumption of contiguity.
 */
template< typename ValueType ,
          class DeviceDst , class MapDst , bool ContigDst ,
          class DeviceSrc , class MapSrc , bool ContigSrc >
class MDArrayDeepCopy {
private:
  enum { DifferentDevices = ! Impl::SameType<DeviceDst,DeviceSrc>::value };
  enum { DifferentMaps    = ! Impl::SameType<MapDst,MapSrc>::value };
  enum { OK = Impl::StaticAssert< DifferentDevices && DifferentMaps >::OK };
public:

  typedef MDArrayView<ValueType,DeviceDst,MapDst> dst_type ;
  typedef MDArrayView<ValueType,DeviceSrc,MapSrc> src_type ;
  typedef MDArrayView<ValueType,DeviceDst,MapSrc> dst_srcmap_type ;

  typedef MDArrayDeepCopy< ValueType,
                           DeviceDst,MapSrc,dst_srcmap_type::Contiguous,
                           DeviceSrc,MapSrc,src_type::Contiguous >
    relocate_operator ;

  typedef MDArrayDeepCopy< ValueType,
                           DeviceDst,MapDst, dst_type::Contiguous,
                           DeviceDst,MapSrc, dst_srcmap_type::Contiguous >
    remap_operator ;

  // Both the devices and the maps are different.
  // Copy to a temporary on the destination with the source map
  // and then remap the temporary to the final array.
  static
  void run( const dst_type & dst , const src_type & src )
  {
    enum { MAX_RANK = 8 };

    size_t dims[ MAX_RANK ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    Impl::mdarray_require_equal_dimension( dst , src );

    dst.dimensions( dims );

    dst_srcmap_type tmp_dst = create_labeled_mdarray<dst_srcmap_type>(
                                "temporary" ,
                                dims[0] , dims[1] , dims[2] , dims[3] ,
                                dims[4] , dims[5] , dims[6] , dims[7] );

    relocate_operator::run( tmp_dst , src );
    remap_operator   ::run( dst ,     tmp_dst );
  }
};

//----------------------------------------------------------------------------
/** \brief  Deep copy with same DeviceType, same Map, and contiguous */

template< typename ValueType , class DeviceType , class MapOpt >
class MDArrayDeepCopy< ValueType , DeviceType , MapOpt , true ,
                                   DeviceType , MapOpt , true >
{
public:
  typedef typename DeviceType::size_type size_type ;

  typedef MDArrayView< ValueType , DeviceType , MapOpt > array_type ;

  typedef Impl::CopyFunctor< ValueType , DeviceType > functor_type ;

  static void run( const array_type & dst , const array_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    parallel_for( dst.size() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );

  }
};

/** \brief  Deep copy with same DeviceType and different Maps */

template< typename ValueType ,
          class DeviceType , class MapDst , bool ContigDst ,
                             class MapSrc , bool ContigSrc >
class MDArrayDeepCopy< ValueType , DeviceType , MapDst , ContigDst ,
                                   DeviceType , MapSrc , ContigSrc >
{
public:
  typedef typename DeviceType::size_type size_type ;

  typedef MDArrayView< ValueType , DeviceType , MapSrc > src_type ;
  typedef MDArrayView< ValueType , DeviceType , MapDst > dst_type ;

  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 8 > deep8 ;
  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 7 > deep7 ;
  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 6 > deep6 ;
  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 5 > deep5 ;
  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 4 > deep4 ;
  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 3 > deep3 ;
  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 2 > deep2 ;
  typedef Impl::MDArrayDeepCopyFunctor< ValueType , DeviceType , MapDst , MapSrc , 1 > deep1 ;

  static
  void run( const dst_type & dst , const src_type & src )
  {
    Impl::mdarray_require_equal_dimension( dst , src );

    const size_t n = dst.size();

    switch ( dst.rank() ) {
    case 8 : parallel_for( n , deep8( dst , src ) ); break ;
    case 7 : parallel_for( n , deep7( dst , src ) ); break ;
    case 6 : parallel_for( n , deep6( dst , src ) ); break ;
    case 5 : parallel_for( n , deep5( dst , src ) ); break ;
    case 4 : parallel_for( n , deep4( dst , src ) ); break ;
    case 3 : parallel_for( n , deep3( dst , src ) ); break ;
    case 2 : parallel_for( n , deep2( dst , src ) ); break ;
    case 1 : parallel_for( n , deep1( dst , src ) ); break ;
    }
  }
};

} // namespace Kokkos

#endif /* KOKKOS_MDARRAYDEEPCOPY_HPP */


