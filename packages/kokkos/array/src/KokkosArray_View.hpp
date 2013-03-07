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

#ifndef KOKKOSARRAY_VIEW_HPP
#define KOKKOSARRAY_VIEW_HPP

#include <string>
#include <KokkosArray_Macros.hpp>

#include <KokkosArray_MemoryManagement.hpp>

#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_Shape.hpp>
#include <impl/KokkosArray_AnalyzeShape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstViewSpecialize , class SrcViewSpecialize = void , class Enable = void >
struct ViewAssignment {};

/** \brief  View specialization mapping of view traits to a specialization tag */
template< class ViewTraits ,
          class ArrayLayout  = typename ViewTraits::array_layout ,
          unsigned Rank      =          ViewTraits::rank ,
          class MemorySpace  = typename ViewTraits::memory_space ,
          class MemoryTraits = typename ViewTraits::memory_traits ,
          class Enable       = void >
struct ViewSpecialize { typedef void type ; };

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType ,
          class LayoutType ,
          class DeviceType ,
          class MemoryTraits >
class ViewTraits ;

/** \brief  ViewTraits
 *
 *  Template argument possibilities:
 *    View< T , Device , Device , void >
 *    View< T , Device , MemoryTraits , void >
 *    View< T , Layout , Device , void >
 *    View< T , Layout , Device , MemoryTraits >
 */

template< class DataType ,
          class LayoutType ,
          class DeviceType ,
          class MemoryTraits >
class ViewTraits {
private:

  typedef Impl::AnalyzeShape<DataType> analysis ;

public:

  typedef DataType                            data_type ;
  typedef LayoutType                          layout_type ;
  typedef DeviceType                          device_type ;
  typedef MemoryTraits                        memory_traits ;

  typedef typename analysis::non_const_type   non_const_data_type ;
  typedef typename analysis::const_type       const_data_type ;
  typedef typename analysis::scalar_type      scalar_type ;
  typedef typename analysis::value_type       value_type ;
  typedef typename analysis::const_value_type const_value_type ;
  typedef typename analysis::shape            shape_type ;
  typedef typename layout_type::array_layout  array_layout ;
  typedef typename device_type::memory_space  memory_space ;
  typedef typename device_type::size_type     size_type ;

  enum { rank         = shape_type::rank };
  enum { rank_dynamic = shape_type::rank_dynamic };
  enum { is_managed   = Impl::is_same< memory_traits , MemoryManaged >::value };
};

/** \brief  Memory traits as third argument, void as fourth argument. */

template< class DataType , class DeviceType , class MemoryTraits >
class ViewTraits< DataType , DeviceType , MemoryTraits , 
  typename Impl::enable_if_type< typename MemoryTraits::memory_traits >::type >
{
private:

  typedef Impl::AnalyzeShape<DataType>       analysis ;

public:

  typedef DataType                            data_type ;
  typedef DeviceType                          layout_type ;
  typedef DeviceType                          device_type ;
  typedef MemoryTraits                        memory_traits ;

  typedef typename analysis::non_const_type   non_const_data_type ;
  typedef typename analysis::const_type       const_data_type ;
  typedef typename analysis::scalar_type      scalar_type ;
  typedef typename analysis::value_type       value_type ;
  typedef typename analysis::const_value_type const_value_type ;
  typedef typename analysis::shape            shape_type ;
  typedef typename layout_type::array_layout  array_layout ;
  typedef typename device_type::memory_space  memory_space ;
  typedef typename device_type::size_type     size_type ;

  enum { rank         = shape_type::rank };
  enum { rank_dynamic = shape_type::rank_dynamic };
  enum { is_managed   = Impl::is_same< memory_traits , MemoryManaged >::value };
};

/** \brief  Device as third argument, void as fourth argument. */

template< class DataType , class LayoutType , class DeviceType >
class ViewTraits< DataType , LayoutType , DeviceType ,
  typename Impl::enable_if_type< typename DeviceType::device_type >::type >
{
private:

  typedef Impl::AnalyzeShape<DataType>       analysis ;

public:

  typedef DataType                            data_type ;
  typedef LayoutType                          layout_type ;
  typedef DeviceType                          device_type ;
  typedef MemoryManaged                       memory_traits ;

  typedef typename analysis::non_const_type   non_const_data_type ;
  typedef typename analysis::const_type       const_data_type ;
  typedef typename analysis::scalar_type      scalar_type ;
  typedef typename analysis::value_type       value_type ;
  typedef typename analysis::const_value_type const_value_type ;
  typedef typename analysis::shape            shape_type ;
  typedef typename layout_type::array_layout  array_layout ;
  typedef typename device_type::memory_space  memory_space ;
  typedef typename device_type::size_type     size_type ;

  enum { rank         = shape_type::rank };
  enum { rank_dynamic = shape_type::rank_dynamic };
  enum { is_managed   = Impl::is_same< memory_traits , MemoryManaged >::value };
};

} // namespace KokkosArray

//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType ,
          class LayoutType ,
          class DeviceType = LayoutType ,
          class MemoryTraits = void ,
          class Specialize   =
            typename Impl::ViewSpecialize< ViewTraits<DataType,LayoutType,DeviceType,MemoryTraits> >::type >
class View ;

//----------------------------------------------------------------------------

template< class LT , class LL , class LD , class LM , class LS ,
          class RT , class RL , class RD , class RM , class RS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator == ( const View<LT,LL,LD,LM,LS> & lhs ,
                   const View<RT,RL,RD,RM,RS> & rhs )
{
  // Same data, layout, dimensions
  typedef ViewTraits<LT,LL,LD,LM> lhs_traits ;
  typedef ViewTraits<RT,RL,RD,RM> rhs_traits ;

  const bool equal_traits =
    Impl::is_same< typename lhs_traits::const_data_type ,
                   typename rhs_traits::const_data_type >::value &&
    Impl::is_same< typename lhs_traits::array_layout ,
                   typename rhs_traits::array_layout >::value &&
    Impl::is_same< typename lhs_traits::memory_space ,
                   typename rhs_traits::memory_space >::value ;

  return equal_traits &&
         lhs.ptr_on_device() == rhs.ptr_on_device() &&
         lhs.dimension_0() == rhs.dimension_0() &&
         lhs.dimension_1() == rhs.dimension_1() &&
         lhs.dimension_2() == rhs.dimension_2() &&
         lhs.dimension_3() == rhs.dimension_3() &&
         lhs.dimension_4() == rhs.dimension_4() &&
         lhs.dimension_5() == rhs.dimension_5() &&
         lhs.dimension_6() == rhs.dimension_6() &&
         lhs.dimension_7() == rhs.dimension_7();
}

template< class LT , class LL , class LD , class LM , class LS ,
          class RT , class RL , class RD , class RM , class RS >
KOKKOSARRAY_INLINE_FUNCTION
bool operator != ( const View<LT,LL,LD,LM,LS> & lhs ,
                   const View<RT,RL,RD,RM,RS> & rhs )
{
  return ! operator==( lhs , rhs );
}

//----------------------------------------------------------------------------

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst , src , arg0 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst, src, arg0, arg1 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst, src, arg0, arg1, arg2 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst, src, arg0, arg1, arg2, arg3 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst, src, arg0, arg1, arg2, arg3, arg4 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst, src, arg0, arg1, arg2, arg3, arg4, arg5 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 ,
         const ArgType6 & arg6 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst, src, arg0, arg1, arg2, arg3, arg4, arg5, arg6 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 , class ArgType7 >
KOKKOSARRAY_INLINE_FUNCTION
View< typename DstViewType::data_type ,
      typename DstViewType::layout_type ,
      typename DstViewType::device_type ,
      MemoryUnmanaged >
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 ,
         const ArgType6 & arg6 ,
         const ArgType7 & arg7 )
{
  typedef ViewTraits< typename DstViewType::data_type ,
                      typename DstViewType::layout_type ,
                      typename DstViewType::device_type ,
                      MemoryUnmanaged > dst_traits ;

  typedef typename Impl::ViewSpecialize<dst_traits>::type dst_spec ;

  typedef View< typename DstViewType::data_type ,
                typename DstViewType::layout_type ,
                typename DstViewType::device_type ,
                MemoryUnmanaged ,
                dst_spec > dst_type ;

  typedef View<T,L,D,M,S> src_type ;

  dst_type dst ;

  Impl::ViewAssignment<dst_spec,S>( dst, src, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7 );

  return dst ;
}

//----------------------------------------------------------------------------


template< class T , class L , class D , class M , class S >
typename View<T,L,D,M,S>::HostMirror
create_mirror_view( const View<T,L,D,M,S> & view ,
                    typename Impl::enable_if<
                      Impl::is_same< typename ViewTraits<T,L,D,M>::memory_space ,
                                     HostSpace >::value
                    >::type * = 0 )
{
  return view ;
}

template< class T , class L , class D , class M , class S >
typename View<T,L,D,M,S>::HostMirror
create_mirror_view( const View<T,L,D,M,S> & view ,
                    typename Impl::enable_if<
                      ! Impl::is_same< typename ViewTraits<T,L,D,M>::memory_space ,
                                       HostSpace >::value
                    >::type * = 0 )
{
  typedef typename View<T,L,D,M>::HostMirror host_view ;
  host_view tmp ;
  Impl::ViewAssignment< S >( tmp , view );
  return tmp ;
}

template< class T , class L , class D , class M , class S >
typename View<T,L,D,M,S>::HostMirror
create_mirror( const View<T,L,D,M,S> & view )
{
#if KOKKOSARRAY_MIRROR_VIEW_OPTIMIZE
  return create_mirror_view( view );
#else
  typedef typename View<T,L,D,M,S>::HostMirror host_view ;
  host_view tmp ;
  Impl::ViewAssignment< S >( tmp , view );
  return tmp ;
#endif
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/KokkosArray_View_factory.hpp>

#include <impl/KokkosArray_ViewAssignment.hpp>

#include <impl/KokkosArray_ViewScalar.hpp>
#include <impl/KokkosArray_ViewVector.hpp>

#include <impl/KokkosArray_ViewLeft.hpp>
#include <impl/KokkosArray_ViewRight.hpp>
#include <impl/KokkosArray_ViewTileLeft.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

