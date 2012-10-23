/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_IMPL_VIEW_FACTORY_HPP
#define KOKKOSARRAY_IMPL_VIEW_FACTORY_HPP

#include <iostream>

#include <KokkosArray_HostSpace.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class ExecutionSpace , class MemoryManagement ,
          class ScalarType , class ShapeType , class MemorySpace >
struct ViewManagement 
{
  KOKKOSARRAY_INLINE_FUNCTION static void increment( ScalarType * ) {}

  KOKKOSARRAY_INLINE_FUNCTION static void decrement( ScalarType * ) {}

  inline static 
  ScalarType * allocate( const std::string & , const ShapeType & )
  { return 0 ; }
};

template< class ScalarType , class ShapeType , class MemorySpace >
struct ViewManagement< HostSpace , MemoryManaged ,
                       ScalarType , ShapeType , MemorySpace >
{
  KOKKOSARRAY_INLINE_FUNCTION static void increment( ScalarType * p )
  { MemorySpace::increment( p ); }

  KOKKOSARRAY_INLINE_FUNCTION static void decrement( ScalarType * p )
  { MemorySpace::decrement( p ); }

  inline static 
  ScalarType * allocate( const std::string & label ,
                         const ShapeType & shape )
  {
    return (ScalarType *)
      MemorySpace::allocate( label ,
                             typeid(ScalarType) ,
                             sizeof(ScalarType) ,
                             Impl::ShapeMap<ShapeType>
                                 ::allocation_count( shape ) );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class > struct ViewCreateMirror ;

template< class DataType , class LayoutType , class DeviceType , class ManageType >
struct ViewCreateMirror< View< DataType , LayoutType , DeviceType , ManageType > >
{
  typedef View< DataType , LayoutType , DeviceType , ManageType > output_type ;

  inline static
  output_type create( const output_type & src ) { return src ; }

  template< class DeviceSrc , class ManageSrc >
  inline static
  output_type create( const View< DataType , LayoutType , DeviceSrc , ManageSrc > & src )
  {
    return output_type( "mirror" , src.shape() );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType ,
          class LayoutType ,
          class DeviceType ,
          class ManageType >
typename View< DataType , LayoutType , DeviceType , ManageType >::HostMirror
create_mirror_view(
  const View<DataType,LayoutType,DeviceType,ManageType > & input )
{
  typedef View< DataType , LayoutType , DeviceType > input_type ;
  typedef typename input_type::HostMirror            output_type ;

  return Impl::ViewCreateMirror< output_type >::create( input );
}

template< class DataType , class LayoutType , class DeviceType >
typename View< DataType , LayoutType , DeviceType >::HostMirror
create_mirror( const View<DataType,LayoutType,DeviceType> & input )
{
  typedef View< DataType , LayoutType , DeviceType > input_type ;
  typedef typename input_type::HostMirror            output_type ;

#if KOKKOSARRAY_MIRROR_VIEW_OPTIMIZE
  return Impl::ViewCreateMirror< output_type >::create( input );
#else
  return output_type( "mirror" , input.shape() );
#endif
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Deep copy compatible arrays */

namespace KokkosArray {
namespace Impl {

template< typename X , typename Y >
struct ViewDeepCopyValueCompatible : public false_type {};

template< typename X >
struct ViewDeepCopyValueCompatible< const X , const X > : public false_type {};

template< typename X >
struct ViewDeepCopyValueCompatible< X , const X > : public true_type {};

template< typename X >
struct ViewDeepCopyValueCompatible< X , X > : public true_type {};


template< class DstView ,
          class SrcView ,
          class SameValue =
            typename ViewDeepCopyValueCompatible<
              typename DstView::scalar_type ,
              typename SrcView::scalar_type >::type ,
          class SameLayout = typename
            is_same< typename DstView::array_layout ,
                     typename SrcView::array_layout >::type ,
          class SameRank = typename
             bool_< DstView::Rank == SrcView::Rank >::type >
struct ViewDeepCopy ;

// Deep copy compatible views:

template< class DataDst , class LayoutDst , class DeviceDst , class ManageDst ,
          class DataSrc , class LayoutSrc , class DeviceSrc , class ManageSrc >
struct ViewDeepCopy< View< DataDst , LayoutDst , DeviceDst , ManageDst > ,
                     View< DataSrc , LayoutSrc , DeviceSrc , ManageSrc > ,
                     true_type /* Same scalar_type   */ ,
                     true_type /* Same array_layout */ ,
                     true_type /* Same rank */ >
{
  typedef View< DataDst , LayoutDst , DeviceDst , ManageDst >  dst_type ;
  typedef View< DataSrc , LayoutSrc , DeviceSrc , ManageSrc >  src_type ;

  inline static
  void apply( const dst_type & dst , const src_type & src )
  {
    typedef typename dst_type::shape_type shape_type ;

    if ( dst != src ) {

      assert_shapes_are_equal( dst.shape() , src.shape() );

      const size_t n =
        sizeof(typename dst_type::scalar_type) *
        ShapeMap< shape_type >::allocation_count( dst.shape() );

        KokkosArray::DeepCopy< typename DeviceDst::memory_space ,
                               typename DeviceSrc::memory_space >(
        dst.ptr_on_device() ,
        src.ptr_on_device() ,
        n );
    }
  }
};

} // namespace Impl

// Deep copy arbitrary arrays:

template< class DataDst , class LayoutDst , class DeviceDst , class ManageDst ,
          class DataSrc , class LayoutSrc , class DeviceSrc , class ManageSrc >
inline
void deep_copy( const View< DataDst, LayoutDst, DeviceDst, ManageDst> & dst ,
                const View< DataSrc, LayoutSrc, DeviceSrc, ManageSrc> & src )
{
  typedef View< DataDst , LayoutDst , DeviceDst , ManageDst > dst_type ;
  typedef View< DataSrc , LayoutSrc , DeviceSrc , ManageSrc > src_type ;

  Impl::ViewDeepCopy<dst_type,src_type>::apply( dst , src );
}

//----------------------------------------------------------------------------

} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

/** \brief  Subview must have compatible pointer types and same memory space */

template< class DstType , class DstMemory ,
          class SrcType , class SrcMemory >
struct SubviewAssignable : public false_type {};

template< class Type , class Memory >
struct SubviewAssignable< Type , Memory , Type , Memory >
  : public true_type {};

template< class Type , class Memory >
struct SubviewAssignable< const Type , Memory , Type , Memory >
  : public true_type {};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_IMPL_VIEW_FACTORY_HPP */

