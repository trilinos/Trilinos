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

#ifndef KOKKOSARRAY_VIEWASSIGNMENT_HPP
#define KOKKOSARRAY_VIEWASSIGNMENT_HPP

#include <stdio.h>

#include <typeinfo>
#include <utility>
#include <KokkosArray_Macros.hpp>

#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class MemorySpace , class MemoryTraits , class = KokkosArray::ExecutionSpace >
struct ViewTracking {
  KOKKOSARRAY_INLINE_FUNCTION static void increment( const void * ) {}
  KOKKOSARRAY_INLINE_FUNCTION static void decrement( const void * ) {}
};

template< class MemorySpace >
struct ViewTracking< MemorySpace , MemoryManaged , HostSpace >
{
  KOKKOSARRAY_INLINE_FUNCTION static void increment( const void * ptr )
    { MemorySpace::increment( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION static void decrement( const void * ptr )
    { MemorySpace::decrement( ptr ); }
};

template< class DstViewType >
struct ViewAssignment< DstViewType , void , unsigned_<0> >
{
  typedef ViewTracking< typename DstViewType::memory_space ,
                        typename DstViewType::memory_traits > tracking ;

  KOKKOSARRAY_INLINE_FUNCTION
  static void increment( const void * ptr ) { tracking::increment( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  static void decrement( const void * ptr ) { tracking::decrement( ptr ); }

  KOKKOSARRAY_INLINE_FUNCTION
  explicit
  ViewAssignment( DstViewType & dst )
  {
    decrement( dst.m_ptr_on_device );
    dst.m_ptr_on_device = 0 ;
  }

  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const DstViewType & src )
  {
    decrement( dst.m_ptr_on_device );

    dst.m_shape         = src.m_shape ;
    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    increment( dst.m_ptr_on_device );
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstShape , class SrcShape ,
          unsigned DstRankDynamic   = DstShape::rank_dynamic ,
          bool     DstRankDynamicOK = DstShape::rank_dynamic >= SrcShape::rank_dynamic >
struct ShapeCompatible { enum { value = false }; };

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 8 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 7 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N7 == SrcShape::N7 };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 6 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N6 == SrcShape::N6 &&
                 DstShape::N7 == SrcShape::N7 };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 5 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N5 == SrcShape::N5 &&
                 DstShape::N6 == SrcShape::N6 &&
                 DstShape::N7 == SrcShape::N7 };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 4 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N4 == SrcShape::N4 &&
                 DstShape::N5 == SrcShape::N5 &&
                 DstShape::N6 == SrcShape::N6 &&
                 DstShape::N7 == SrcShape::N7 };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 3 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N3 == SrcShape::N3 &&
                 DstShape::N4 == SrcShape::N4 &&
                 DstShape::N5 == SrcShape::N5 &&
                 DstShape::N6 == SrcShape::N6 &&
                 DstShape::N7 == SrcShape::N7 };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 2 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N2 == SrcShape::N2 &&
                 DstShape::N3 == SrcShape::N3 &&
                 DstShape::N4 == SrcShape::N4 &&
                 DstShape::N5 == SrcShape::N5 &&
                 DstShape::N6 == SrcShape::N6 &&
                 DstShape::N7 == SrcShape::N7 };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 1 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N1 == SrcShape::N1 &&
                 DstShape::N2 == SrcShape::N2 &&
                 DstShape::N3 == SrcShape::N3 &&
                 DstShape::N4 == SrcShape::N4 &&
                 DstShape::N5 == SrcShape::N5 &&
                 DstShape::N6 == SrcShape::N6 &&
                 DstShape::N7 == SrcShape::N7 };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 0 , true >
{
  enum { value = DstShape::scalar_size == SrcShape::scalar_size &&
                 DstShape::N0 == SrcShape::N0 &&
                 DstShape::N1 == SrcShape::N1 &&
                 DstShape::N2 == SrcShape::N2 &&
                 DstShape::N3 == SrcShape::N3 &&
                 DstShape::N4 == SrcShape::N4 &&
                 DstShape::N5 == SrcShape::N5 &&
                 DstShape::N6 == SrcShape::N6 &&
                 DstShape::N7 == SrcShape::N7 };
};

//----------------------------------------------------------------------------

template< class DstViewTraits , class SrcViewTraits >
struct ViewAssignment_Compatible {
public:

  enum { compatible_value =
    is_same< typename DstViewTraits::memory_space ,
             typename SrcViewTraits::memory_space >::value &&
    ( is_same< typename DstViewTraits::value_type ,
               typename SrcViewTraits::value_type >::value ||
      is_same< typename DstViewTraits::value_type ,
               typename SrcViewTraits::const_value_type >::value )
  };

  enum { compatible_layout = is_same< typename DstViewTraits::array_layout ,
                                      typename SrcViewTraits::array_layout >::value };

  enum { compatible_shape = ShapeCompatible< typename DstViewTraits::shape_type ,
                                             typename SrcViewTraits::shape_type >::value };

  enum { compatible = compatible_value && compatible_layout && compatible_shape };
};

//----------------------------------------------------------------------------
// 1) Fully compatible DstViewType and SrcViewType
// 2) ArgCount = 0

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ViewAssignment_Compatible< typename DstViewType::view_traits ,
                               typename SrcViewType::view_traits >::compatible
  ) , unsigned_<0> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src )
  {
    typedef typename DstViewType::shape_type shape_type ;

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape ,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 1 )
    &&
    ( SrcViewType::Rank == 1 )
  ) , unsigned_<1> >::type >
{
  template< typename iType >
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const std::pair<iType,iType> & range )
  {
    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    if ( range.first < range.second ) {

      assert_shape_bounds( src.shape() , range.first );
      assert_shape_bounds( src.shape() , range.second - 1 );

      dst.m_shape.N0 = range.second - range.first ;
      dst.m_stride   = 0 ;
      dst.m_ptr_on_device = src.m_ptr_on_device + range.first ;

      ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
    }
    else {
      dst.m_shape.N0 = 0 ;
      dst.m_stride   = 0 ;
      dst.m_ptr_on_device = 0 ;
    }
  }
};

//----------------------------------------------------------------------------
// Vector from multivector:

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 1 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutLeft >::value )
    &&
    ( SrcViewType::Rank == 2 )
  ) , unsigned_<1> >::type >
{
  inline
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i1 )
  {
    assert_shape_bounds( src.shape() , 0 , i1 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0 = src.m_shape.N0 ;
    dst.m_stride   = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device + src.m_stride * i1 ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

#if 0
//----------------------------------------------------------------------------
// ArgCount=2 , dst<Rank=1> : src<Rank=2,LayoutLeft>  : ( range , index )
// ArgCount=2 , dst<Rank=1> : src<Rank=2,LayoutRight> : ( index , range )
//----------------------------------------------------------------------------

// Assignment of rank-zero array

template< class DstViewType >
struct ViewAssignment< DstViewType , 0 >
{
  template< const SrcViewType >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  typename enable_if<
                    ViewAssignment_Compatible< typename DstViewType::view_traits ,
                                               typename SrcViewType::view_traits >::compatible_value
                  >::type * = 0 )
  {
    typedef typename DstDeviceType::memory_space dst_memory_space ;

    if ( DstManagement::managed ) dst_memory_space::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    if ( DstManagement::managed ) dst_memory_space::increment( dst.m_ptr_on_device );
  }

  template< const SrcViewType , typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const iType0 & i0 ,
                  typename enable_if<
                    ViewAssignment_Compatible< typename DstViewType::view_traits ,
                                               typename SrcViewType::view_traits >::compatible_value
                  >::type * = 0 )
  {
    typedef typename DstDeviceType::memory_space dst_memory_space ;

    if ( DstManagement::managed ) dst_memory_space::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = & src(i0);

    if ( DstManagement::managed ) dst_memory_space::increment( dst.m_ptr_on_device );
  }

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const iType0 & i0 , const iType1 & i1 )
  {
    typedef typename DstDeviceType::memory_space  dst_memory_space ;

    if ( DstManagement::managed ) dst_memory_space::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = & src(i0,i1);

    if ( DstManagement::managed ) dst_memory_space::increment( dst.m_ptr_on_device );
  }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const iType0 & i0 , const iType1 & i1 , const iType2 & i2 )
  {
    typedef typename DstDeviceType::memory_space  dst_memory_space ;

    if ( DstManagement::managed ) dst_memory_space::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = & src(i0,i1,i2);

    if ( DstManagement::managed ) dst_memory_space::increment( dst.m_ptr_on_device );
  }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 )
  {
    typedef typename DstDeviceType::memory_space  dst_memory_space ;

    if ( DstManagement::managed ) dst_memory_space::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = & src(i0,i1,i2,i3);

    if ( DstManagement::managed ) dst_memory_space::increment( dst.m_ptr_on_device );
  }




  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                  const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 )
  {
    typedef typename DstDeviceType::memory_space  dst_memory_space ;

    if ( DstManagement::managed ) dst_memory_space::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = & src(i0,i1,i2,i3,i4,i5,i6,i7);

    if ( DstManagement::managed ) dst_memory_space::increment( dst.m_ptr_on_device );
  }
};

#endif

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_VIEWASSIGNMENT_HPP */

