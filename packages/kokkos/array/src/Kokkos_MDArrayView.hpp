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

#ifndef KOKKOS_MDARRAYVIEW_HPP
#define KOKKOS_MDARRAYVIEW_HPP

#include <cstddef>
#include <string>
#include <Kokkos_ArrayForwardDeclarations.hpp>
#include <impl/Kokkos_ArrayBounds.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

namespace Kokkos {

template< typename ValueType , class DeviceType , class MapOption >
MDArrayView< ValueType , DeviceType , MapOption >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 );

namespace Impl {

template< typename ValueType ,
          class DeviceDst , class MapDst , bool ContigDst ,
          class DeviceSrc , class MapSrc , bool ContigSrc >
class MDArrayDeepCopy ;

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Multidimensional array allocated and mapped
 *          onto a compute device.
 *
 *  The array is a simple rank-N container of simple scalar values
 *  where 1 <= N <= 8.
 *
 *  The first rank is the parallel work index.
 *
 *  No assumptions should be made as to the mapping, contiguity, or strides
 *  of the storage of these arrays.  The mapping will vary according to the
 *  underlying device.  The usage model is for algorithms to be parameterized
 *  with respect to the type of the mapped array and thus achieve portability
 *  across compute devices.
 *
 *  Several APIs for MDArray creation functions are available.
 *  The "labeled" group creates MDArray's with string labels that
 *  will appear in error messages.
 *
 *  create_labeled_mdarray< ValueType , DeviceType , MapOption >( label , nP , ... );
 *  create_labeled_mdarray< ValueType , DeviceType >( label , nP , ... );
 *  create_labeled_mdarray< MDArrayView<...> >( label , nP , ... );
 *
 *  The "unlabeled" group creates MDArray's with NULL string labels.
 *
 *  create_mdarray< ValueType , DeviceType , MapOption >( nP , ... );
 *  create_mdarray< ValueType , DeviceType >( nP , ... );
 *  create_mdarray< MDArrayView<...> >( nP , ... );
 */

template< typename ValueType , class DeviceType , class MapOption >
class MDArrayView {
public:
  typedef ValueType                       value_type ;
  typedef DeviceType                      device_type ;
  typedef MapOption                       map_option ;
  typedef typename DeviceType::size_type  size_type ;

  /*------------------------------------------------------------------*/
  /** \brief  True if the array type has contigous memory
   *          If contigous then can get a pointer to the memory.
   */
  enum { Contiguous = false };

  /** \brief  Query rank of the array */
  size_type rank() const ;

  /** \brief  Query dimension of the given ordinate of the array */
  template < typename iType >
  size_type dimension( const iType & rank_ordinate ) const ;

  /** \brief  Query all dimensions */
  template< typename iType >
  void dimensions( iType * const dims ) const ;

  /** \brief  Query total number of members */
  size_type size() const ;

  /*------------------------------------------------------------------*/
  /** \brief  Query value of a rank 8 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 , const iType7 & i7 ) const ;

  /** \brief  Query value of a rank 7 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ,
                           const iType6 & i6 ) const ;

  /** \brief  Query value of a rank 6 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 , const iType5 & i5 ) const ;

  /** \brief  Query value of a rank 5 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ,
                           const iType4 & i4 ) const ;

  /** \brief  Query value of a rank 4 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 , typename iType3 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 , const iType3 & i3 ) const ;

  /** \brief  Query value of a rank 3 array */
  template< typename iTypeP , typename iType1 ,
            typename iType2 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ,
                           const iType2 & i2 ) const ;

  /** \brief  Query value of a rank 2 array */
  template< typename iTypeP , typename iType1 >
  value_type & operator()( const iTypeP & iP , const iType1 & i1 ) const ;

  /** \brief  Query value of a rank 1 array */
  template< typename iTypeP >
  value_type & operator()( const iTypeP & iP ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  MDArrayView();

  /** \brief  Construct another view of the 'rhs' array */
  inline
  MDArrayView( const MDArrayView & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  MDArrayView & operator = ( const MDArrayView & rhs );
  
  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  ~MDArrayView();

  /*------------------------------------------------------------------*/
  /** \brief  Query if non-NULL view */
  operator bool () const ;

  bool operator == ( const MDArrayView & ) const ;

  bool operator != ( const MDArrayView & ) const ;

private:

  MDArrayView( const std::string & label ,
               size_t nP , size_t n1 , size_t n2 , size_t n3 ,
               size_t n4 , size_t n5 , size_t n6 , size_t n7 );

  template< typename V , class D , class M >
  friend
  MDArrayView< V , D , M >
  create_labeled_mdarray( const std::string & label ,
                          size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                          size_t n4 , size_t n5 , size_t n6 , size_t n7 );

  
  template< typename V , class DeviceDst , class MapDst , bool ContigDst ,
                         class DeviceSrc , class MapSrc , bool ContigSrc >
  friend
  class Impl::MDArrayDeepCopy ;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  This is THE creation function.
 *          All other versions call this function.
 */
template< typename ValueType , class DeviceType , class MapType >
inline
MDArrayView< ValueType , DeviceType , MapType >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 , size_t n2 , size_t n3 ,
                        size_t n4 , size_t n5 , size_t n6 , size_t n7 )
{
  return MDArrayView< ValueType , DeviceType , MapType >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
MDArrayView< ValueType , DeviceType , typename DeviceType::default_map >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 )
{
  return create_labeled_mdarray< ValueType , DeviceType , typename DeviceType::default_map >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename MDArrayType >
inline
MDArrayView< typename MDArrayType::value_type ,
             typename MDArrayType::device_type ,
             typename MDArrayType::map_option >
create_labeled_mdarray( const std::string & label ,
                        size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                                    size_t n3 = 0 , size_t n4 = 0 ,
                                    size_t n5 = 0 , size_t n6 = 0 ,
                                    size_t n7 = 0 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type ,
                                 typename MDArrayType::map_option >
           ( label , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType , class MapType >
inline
MDArrayView< ValueType , DeviceType , MapType >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< ValueType , DeviceType , MapType >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
MDArrayView< ValueType , DeviceType , typename DeviceType::default_map >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< ValueType , DeviceType , typename DeviceType::default_map >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename MDArrayType >
inline
MDArrayView< typename MDArrayType::value_type ,
             typename MDArrayType::device_type ,
             typename MDArrayType::map_option >
create_mdarray( size_t nP , size_t n1 = 0 , size_t n2 = 0 ,
                            size_t n3 = 0 , size_t n4 = 0 ,
                            size_t n5 = 0 , size_t n6 = 0 ,
                            size_t n7 = 0 )
{
  return create_labeled_mdarray< typename MDArrayType::value_type ,
                                 typename MDArrayType::device_type ,
                                 typename MDArrayType::map_option >
           ( std::string() , nP , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
}

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceDst , class MapDst ,
                               class DeviceSrc , class MapSrc >
inline
void deep_copy( const MDArrayView<ValueType,DeviceDst,MapDst> & dst ,
                const MDArrayView<ValueType,DeviceSrc,MapSrc> & src )
{
  enum { ContigDst = MDArrayView<ValueType,DeviceDst,MapDst>::Contiguous };
  enum { ContigSrc = MDArrayView<ValueType,DeviceSrc,MapSrc>::Contiguous };

  typedef MDArrayView<ValueType,DeviceDst,MapDst>  dst_type ;
  typedef MDArrayView<ValueType,DeviceSrc,MapSrc>  src_type ;

  Impl::mdarray_require_equal_dimension( dst , src );

  Impl::MDArrayDeepCopy< ValueType,
                         DeviceDst, MapDst, ContigDst,
                         DeviceSrc, MapSrc, ContigSrc >::run( dst , src );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Partial specializations for known devices and default index maps

#if defined( KOKKOS_DEVICE_HOST )
#include <Kokkos_DeviceHost_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <impl/Kokkos_MDArrayView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_TPI )
#include <Kokkos_DeviceTPI_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <impl/Kokkos_MDArrayView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>
#include <DeviceTPI/Kokkos_DeviceTPI_MDArrayView.hpp>
#endif

#if defined( KOKKOS_DEVICE_CUDA )
#include <Kokkos_DeviceCuda_macros.hpp>
#include <impl/Kokkos_MDArrayIndexMapLeft_macros.hpp>
#include <impl/Kokkos_MDArrayView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>
#include <DeviceCuda/Kokkos_DeviceCuda_MDArrayView.hpp>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Deep copy between different devices, different maps,
 *          and no assumption of contiguity.
 */
template< typename ValueType , class DeviceDst , class MapDst , bool ,
                               class DeviceSrc , class MapSrc , bool >
class MDArrayDeepCopy {
private:
  enum { okD = StaticAssert< ! SameType<DeviceDst,DeviceSrc>::value >::value };
  enum { okM = StaticAssert< ! SameType<MapDst,MapSrc>::value >::value };
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
    size_t dims[ MDArrayMaxRank ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

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

} // namespace Impl
} // namespace Kokkos

#endif /* KOKKOS_MDARRAYVIEW_HPP */


