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

#include <typeinfo>
#include <utility>
#include <KokkosArray_Macros.hpp>

#include <impl/KokkosArray_ArrayTraits.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class DstViewType >
struct ViewInitialize { static void apply( const DstViewType & ) {} };

template< class T , class L , class D , class M , class S >
size_t allocation_count( const View<T,L,D,M,S> & view )
{
  return ViewAssignment<S>::allocation_count( view );
}

//----------------------------------------------------------------------------

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

//----------------------------------------------------------------------------

template< class DstViewType >
struct ViewAssignment< DstViewType , void , void >
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
          bool     DstRankDynamicOK = unsigned(DstShape::rank_dynamic) >= unsigned(SrcShape::rank_dynamic) >
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

template< class DstView , class SrcView ,
          class DstValueType  = typename DstView::value_type ,
          class DstValueSpace = typename DstView::memory_space ,
          class SrcValueType  = typename SrcView::value_type ,
          class SrcValueSpace = typename SrcView::memory_space >
struct ValueCompatible ;

template< class DstView , class SrcView , class ValueType , class ValueSpace >
struct ValueCompatible< DstView , SrcView ,
                        ValueType , ValueSpace ,
                        ValueType , ValueSpace >
{
  typedef ValueType type ;
  enum { value = true };
};

template< class DstView , class SrcView , class ValueType , class ValueSpace >
struct ValueCompatible< DstView , SrcView ,
                        const ValueType , ValueSpace ,
                              ValueType , ValueSpace >
{
  typedef ValueType type ;
  enum { value = true };
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

template< class DT, class DL, class DD, class DM, class Spec,
          class ST, class SL, class SD, class SM >
struct ViewAssignment<
  View<DT,DL,DD,DM,Spec> ,
  View<ST,SL,SD,SM,Spec> ,
  typename enable_if<(
    ViewAssignment_Compatible< ViewTraits<DT,DL,DD,DM> ,
                               ViewTraits<ST,SL,SD,SM> >::compatible
  ), unsigned_<0> >::type >
{
  typedef View<DT,DL,DD,DM,Spec> DstViewType ;
  typedef View<ST,SL,SD,SM,Spec> SrcViewType ;

  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src )
  { ViewAssignment<Spec,Spec>( dst , src ); }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_VIEWASSIGNMENT_HPP */

