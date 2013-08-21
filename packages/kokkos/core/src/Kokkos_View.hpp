/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef KOKKOS_VIEW_HPP
#define KOKKOS_VIEW_HPP

#include <string>
#include <Kokkos_Macros.hpp>
#include <Kokkos_HostSpace.hpp>

#include <Kokkos_MemoryTraits.hpp>

#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_ArrayTraits.hpp>
#include <impl/Kokkos_Shape.hpp>
#include <impl/Kokkos_AnalyzeShape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  View specialization mapping of view traits to a specialization tag */
template< typename ScalarType , class ValueType ,
          class ArrayLayout , class uRank , class uRankDynamic ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize ;

} /* namespace Impl */

/** \brief  ViewTraits
 *
 *  Template argument permutations:
 *
 *    View< DataType , Device , void         , void >
 *    View< DataType , Device , MemoryTraits , void >
 *    View< DataType , Device , void         , MemoryTraits >
 *    View< DataType , ArrayLayout , Device  , void >
 *    View< DataType , ArrayLayout , Device  , MemoryTraits >
 */

template< class DataType ,
          class ArrayLayout ,
          class DeviceType ,
          class MemoryTraits >
class ViewTraits {
private:

  typedef Impl::AnalyzeShape<DataType> analysis ;

public:

  //------------------------------------
  // Data type traits:

  typedef DataType                            data_type ;
  typedef typename analysis::const_type       const_data_type ;
  typedef typename analysis::non_const_type   non_const_data_type ;

  //------------------------------------
  // Scalar type traits:

  typedef typename analysis::scalar_type            scalar_type ;
  typedef typename analysis::const_scalar_type      const_scalar_type ;
  typedef typename analysis::non_const_scalar_type  non_const_scalar_type ;

  //------------------------------------
  // Value type traits:

  typedef typename analysis::value_type            value_type ;
  typedef typename analysis::const_value_type      const_value_type ;
  typedef typename analysis::non_const_value_type  non_const_value_type ;

  //------------------------------------
  // Layout and shape traits:

  typedef typename Impl::StaticAssertSame< ArrayLayout , typename ArrayLayout ::array_layout >::type  array_layout ;

  typedef typename analysis::shape   shape_type ;

  enum { rank         = shape_type::rank };
  enum { rank_dynamic = shape_type::rank_dynamic };

  //------------------------------------
  // Device and memory space traits:

  typedef typename Impl::StaticAssertSame< DeviceType   , typename DeviceType  ::device_type   >::type  device_type ;
  typedef typename Impl::StaticAssertSame< MemoryTraits , typename MemoryTraits::memory_traits >::type  memory_traits ;

  typedef typename device_type::memory_space  memory_space ;
  typedef typename device_type::size_type     size_type ;

  enum { is_hostspace = Impl::is_same< memory_space , HostSpace >::value };
  enum { is_managed   = memory_traits::Unmanaged == 0 };

  //------------------------------------
  // Specialization:
  typedef typename
    Impl::ViewSpecialize< scalar_type ,
                          value_type ,
                          array_layout ,
                          Impl::unsigned_<rank> ,
                          Impl::unsigned_<rank_dynamic> ,
                          memory_space ,
                          memory_traits
                        >::type specialize ;
};

/** \brief  Traits for View<DataType,DeviceType,void,void> */

template< class DataType ,
          class DeviceType >
class ViewTraits<DataType,DeviceType,void,void>
  : public ViewTraits< DataType , typename DeviceType::array_layout , DeviceType , MemoryManaged > {};

/** \brief  Traits for View<DataType,DeviceType,void,MemoryTraits> */

template< class DataType ,
          class DeviceType ,
          class MemoryTraits >
class ViewTraits<DataType,DeviceType,void,MemoryTraits>
  : public ViewTraits< DataType , typename DeviceType::array_layout , DeviceType , MemoryTraits > {};

/** \brief  Traits for View<DataType,DeviceType,MemoryTraits,void> */

template< class DataType , class DeviceType , class MemoryTraits >
class ViewTraits< DataType , DeviceType , MemoryTraits ,
  typename Impl::enable_if< (
    Impl::is_same< DeviceType   , typename DeviceType  ::device_type   >::value &&
    Impl::is_same< MemoryTraits , typename MemoryTraits::memory_traits >::value
  ) >::type >
  : public ViewTraits< DataType , typename DeviceType::array_layout , DeviceType , MemoryTraits > {};

/** \brief  Traits for View<DataType,ArrayLayout,DeviceType,void> */

template< class DataType , class ArrayLayout , class DeviceType >
class ViewTraits< DataType , ArrayLayout , DeviceType ,
  typename Impl::enable_if< (
    Impl::is_same< ArrayLayout , typename ArrayLayout::array_layout >::value &&
    Impl::is_same< DeviceType  , typename DeviceType ::device_type  >::value
  ) >::type >
  : public ViewTraits< DataType , ArrayLayout , DeviceType , MemoryManaged > {};

//----------------------------------------------------------------------------

/** \brief  View to array of data.
 *
 *  Options for template arguments:
 *
 *    View< DataType , Device >
 *    View< DataType , Device ,        MemoryTraits >
 *    View< DataType , Device , void , MemoryTraits >
 *
 *    View< DataType , Layout , Device >
 *    View< DataType , Layout , Device , MemoryTraits >
 */

template< class DataType ,
          class Arg1Type ,        /* ArrayLayout or DeviceType */
          class Arg2Type = void , /* DeviceType or MemoryTraits */
          class Arg3Type = void , /* MemoryTraits */
          class Specialize =
            typename ViewTraits<DataType,Arg1Type,Arg2Type,Arg3Type>::specialize >
class View ;

//----------------------------------------------------------------------------

template< class LT , class LL , class LD , class LM , class LS ,
          class RT , class RL , class RD , class RM , class RS >
KOKKOS_INLINE_FUNCTION
typename Impl::enable_if<( Impl::is_same< LS , RS >::value ), bool >::type
operator == ( const View<LT,LL,LD,LM,LS> & lhs ,
              const View<RT,RL,RD,RM,RS> & rhs )
{
  // Same data, layout, dimensions
  typedef ViewTraits<LT,LL,LD,LM> lhs_traits ;
  typedef ViewTraits<RT,RL,RD,RM> rhs_traits ;

  return
    Impl::is_same< typename lhs_traits::const_data_type ,
                   typename rhs_traits::const_data_type >::value &&
    Impl::is_same< typename lhs_traits::array_layout ,
                   typename rhs_traits::array_layout >::value &&
    Impl::is_same< typename lhs_traits::memory_space ,
                   typename rhs_traits::memory_space >::value &&
    Impl::is_same< typename lhs_traits::specialize ,
                   typename rhs_traits::specialize >::value &&
    lhs.ptr_on_device() == rhs.ptr_on_device() &&
    lhs.shape()         == rhs.shape() ;
}

template< class LT , class LL , class LD , class LM , class LS ,
          class RT , class RL , class RD , class RM , class RS >
KOKKOS_INLINE_FUNCTION
bool operator != ( const View<LT,LL,LD,LM,LS> & lhs ,
                   const View<RT,RL,RD,RM,RS> & rhs )
{
  return ! operator==( lhs , rhs );
}

//----------------------------------------------------------------------------

/** \brief  Reallocate a view without copying old data to new data */
template< class T , class L , class D , class M , class S >
inline
void realloc( View<T,L,D,M,S> & v ,
              const typename Impl::enable_if< ViewTraits<T,L,D,M>::is_managed , size_t >::type n0 ,
              const size_t n1 = 0 ,
              const size_t n2 = 0 ,
              const size_t n3 = 0 ,
              const size_t n4 = 0 ,
              const size_t n5 = 0 ,
              const size_t n6 = 0 ,
              const size_t n7 = 0 )
{
  typedef View<T,L,D,M,S> view_type ;
  typedef typename view_type::memory_space memory_space ;

  // Query the current label and reuse it.
  const std::string label = memory_space::query_label( v.ptr_on_device() );

  v = view_type(); // deallocate first, if the only view to memory.
  v = view_type( label, n0, n1, n2, n3, n4, n5, n6, n7 );
}


//----------------------------------------------------------------------------

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class DstViewSpecialize , class SrcViewSpecialize = void , class Enable = void >
struct ViewAssignment ;

/** \brief  Evaluate if LHS = RHS view assignment is allowed. */
template< class ViewLHS , class ViewRHS >
struct ViewAssignable
{
  // Same memory space.
  // Same value type.
  // Compatible 'const' qualifier
  // Cannot assign managed = unmannaged
  enum { assignable_value =
    ( is_same< typename ViewLHS::value_type ,
               typename ViewRHS::value_type >::value
      ||
      is_same< typename ViewLHS::value_type ,
               typename ViewRHS::const_value_type >::value )
    &&
    is_same< typename ViewLHS::memory_space ,
             typename ViewRHS::memory_space >::value
    &&
    ( ! ( ViewLHS::is_managed && ! ViewRHS::is_managed ) )
  };

  enum { assignable_shape =
    // Compatible shape and matching layout:
    ( ShapeCompatible< typename ViewLHS::shape_type ,
                       typename ViewRHS::shape_type >::value
      &&
      is_same< typename ViewLHS::array_layout ,
               typename ViewRHS::array_layout >::value )
    ||
    // Matching layout, same rank, and LHS dynamic rank
    ( is_same< typename ViewLHS::array_layout ,
               typename ViewRHS::array_layout >::value
      &&
      int(ViewLHS::rank) == int(ViewRHS::rank)
      &&
      int(ViewLHS::rank) == int(ViewLHS::rank_dynamic) )
    ||
    // Both rank-0, any shape and layout
    ( int(ViewLHS::rank) == 0 && int(ViewRHS::rank) == 0 )
    ||
    // Both rank-1 and LHS is dynamic rank-1, any shape and layout
    ( int(ViewLHS::rank) == 1 && int(ViewRHS::rank) == 1 &&
      int(ViewLHS::rank_dynamic) == 1 )
    };

  enum { value = assignable_value && assignable_shape };
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class OutputView , unsigned Rank = OutputView::Rank >
struct ViewInit
{
  typedef typename OutputView::device_type device_type ;
  typedef typename OutputView::scalar_type scalar_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;

  explicit ViewInit( const OutputView & arg_out ) : output( arg_out )
    { parallel_for( output.dimension_0() , *this ); }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    const scalar_type default_value = scalar_type();

    for ( size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      new (&output(i0,i1,i2,i3,i4,i5,i6,i7)) scalar_type(default_value) ;
    }}}}}}}
  }
};

template< class OutputView >
struct ViewInit< OutputView , 1 >
{
  typedef typename OutputView::device_type device_type ;
  typedef typename OutputView::value_type  value_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;

  explicit ViewInit( const OutputView & arg_out ) : output( arg_out )
    { parallel_for( output.dimension_0() , *this ); }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    value_type default_value = value_type();
    new (&output(i0)) value_type(default_value) ;
  }
};

template< class OutputView >
struct ViewInit< OutputView , 0 >
{
  typedef typename OutputView::device_type device_type ;
  typedef typename OutputView::value_type  value_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;

  explicit ViewInit( const OutputView & arg_out ) : output( arg_out )
    { parallel_for( 1 , *this ); }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type /*i0*/ ) const
  {
    value_type default_value = value_type();
    new (&(*output)) value_type(default_value) ;
  }
};

template< class Device >
struct ViewInitialize
{
  template< class ViewType >
  inline explicit ViewInitialize( const ViewType & view )
    { ViewInit<ViewType> init( view ); }
};

template< class OutputView , class InputView  , unsigned Rank = OutputView::Rank >
struct ViewRemap
{
  typedef typename OutputView::device_type device_type ;
  typedef typename device_type::size_type  size_type ;

  const OutputView output ;
  const InputView  input ;
  const size_type n0 ;
  const size_type n1 ;
  const size_type n2 ;
  const size_type n3 ;
  const size_type n4 ;
  const size_type n5 ;
  const size_type n6 ;
  const size_type n7 ;

  ViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    , n0( std::min( (size_t)arg_out.dimension_0() , (size_t)arg_in.dimension_0() ) )
    , n1( std::min( (size_t)arg_out.dimension_1() , (size_t)arg_in.dimension_1() ) )
    , n2( std::min( (size_t)arg_out.dimension_2() , (size_t)arg_in.dimension_2() ) )
    , n3( std::min( (size_t)arg_out.dimension_3() , (size_t)arg_in.dimension_3() ) )
    , n4( std::min( (size_t)arg_out.dimension_4() , (size_t)arg_in.dimension_4() ) )
    , n5( std::min( (size_t)arg_out.dimension_5() , (size_t)arg_in.dimension_5() ) )
    , n6( std::min( (size_t)arg_out.dimension_6() , (size_t)arg_in.dimension_6() ) )
    , n7( std::min( (size_t)arg_out.dimension_7() , (size_t)arg_in.dimension_7() ) )
    {
      parallel_for( n0 , *this );
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type i0 ) const
  {
    for ( size_type i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_type i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_type i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_type i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_type i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_type i6 = 0 ; i6 < n6 ; ++i6 ) {
    for ( size_type i7 = 0 ; i7 < n7 ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}
  }
};

template< class OutputView , class InputView  >
struct ViewRemap< OutputView ,  InputView , 0 >
{
  typedef typename OutputView::value_type   value_type ;
  typedef typename OutputView::device_space dst_space ;
  typedef typename InputView ::device_space src_space ;

  ViewRemap( const OutputView & arg_out , const InputView & arg_in )
  {
    DeepCopy< dst_space , src_space >( arg_out.ptr_on_device() ,
                                       arg_in.ptr_on_device() ,
                                       sizeof(value_type) );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of the same specialization, compatible type,
 *          same rank, same layout are handled by that specialization.
 */

template< class DT , class DL , class DD , class DM , class DS ,
          class ST , class SL , class SD , class SM , class SS >
inline
void deep_copy( const View<DT,DL,DD,DM,DS> & dst ,
                const View<ST,SL,SD,SM,SS> & src ,
                typename Impl::enable_if<(
                  Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::scalar_type ,
                                 typename ViewTraits<ST,SL,SD,SM>::non_const_scalar_type >::value
                  &&
                  Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::array_layout ,
                                 typename ViewTraits<ST,SL,SD,SM>::array_layout >::value
                  &&
                  ( unsigned(ViewTraits<DT,DL,DD,DM>::rank) == unsigned(ViewTraits<ST,SL,SD,SM>::rank) )
                )>::type * = 0 )
{ Impl::ViewAssignment<DS,SS>::deep_copy( dst , src ); }


/** \brief Deep copy equal dimension arrays in the host space which
 *         have different layouts or specializations.
 */
template< class DT , class DL , class DD , class DM , class DS ,
          class ST , class SL ,            class SM , class SS >
inline
void deep_copy( const View< DT, DL, DD, DM, DS> & dst ,
                const View< ST, SL, DD, SM, SS> & src ,
                const typename Impl::enable_if<(
                  // Destination is not constant:
                  Impl::is_same< typename ViewTraits<DT,DL,DD,DM>::value_type ,
                                 typename ViewTraits<DT,DL,DD,DM>::non_const_value_type >::value
                  &&
                  // Same rank
                  ( unsigned( ViewTraits<DT,DL,DD,DM>::rank ) ==
                    unsigned( ViewTraits<ST,SL,DD,SM>::rank ) )
                  &&
                  // Different layout or different specialization:
                  ( ( ! Impl::is_same< typename DL::array_layout ,
                                       typename SL::array_layout >::value )
                    ||
                    ( ! Impl::is_same< DS , SS >::value )
                  )
                )>::type * = 0 )
{
  typedef View< DT, DL, DD, DM, DS> dst_type ;
  typedef View< ST, SL, DD, SM, SS> src_type ;

  assert_shapes_equal_dimension( dst.shape() , src.shape() );

  Impl::ViewRemap< dst_type , src_type >( dst , src );
}

//----------------------------------------------------------------------------

/** \brief  Resize a view with copying old data to new data at the corresponding indices. */
template< class T , class L , class D , class M , class S >
inline
void resize( View<T,L,D,M,S> & v ,
             const typename Impl::enable_if< ViewTraits<T,L,D,M>::is_managed , size_t >::type n0 ,
             const size_t n1 = 0 ,
             const size_t n2 = 0 ,
             const size_t n3 = 0 ,
             const size_t n4 = 0 ,
             const size_t n5 = 0 ,
             const size_t n6 = 0 ,
             const size_t n7 = 0 )
{
  typedef View<T,L,D,M,S> view_type ;
  typedef typename view_type::memory_space memory_space ;

  const std::string label = memory_space::query_label( v.ptr_on_device() );

  view_type v_resized( label, n0, n1, n2, n3, n4, n5, n6, n7 );

  Impl::ViewRemap< view_type , view_type >( v_resized , v );

  v = v_resized ;
}

//----------------------------------------------------------------------------

struct ALL { KOKKOS_INLINE_FUNCTION ALL(){} };

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 >
KOKKOS_INLINE_FUNCTION
DstViewType
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 )
{
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst , src , arg0 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 >
KOKKOS_INLINE_FUNCTION
DstViewType
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 )
{
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst, src, arg0, arg1 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 >
KOKKOS_INLINE_FUNCTION
DstViewType
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 )
{
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst, src, arg0, arg1, arg2 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 >
KOKKOS_INLINE_FUNCTION
DstViewType
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 )
{
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst, src, arg0, arg1, arg2, arg3 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 >
KOKKOS_INLINE_FUNCTION
DstViewType
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 )
{
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst, src, arg0, arg1, arg2, arg3, arg4 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 >
KOKKOS_INLINE_FUNCTION
DstViewType
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 )
{
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst, src, arg0, arg1, arg2, arg3, arg4, arg5 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 >
KOKKOS_INLINE_FUNCTION
DstViewType
subview( const View<T,L,D,M,S> & src ,
         const ArgType0 & arg0 ,
         const ArgType1 & arg1 ,
         const ArgType2 & arg2 ,
         const ArgType3 & arg3 ,
         const ArgType4 & arg4 ,
         const ArgType5 & arg5 ,
         const ArgType6 & arg6 )
{
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst, src, arg0, arg1, arg2, arg3, arg4, arg5, arg6 );

  return dst ;
}

template< class DstViewType ,
          class T , class L , class D , class M , class S ,
          class ArgType0 , class ArgType1 , class ArgType2 , class ArgType3 ,
          class ArgType4 , class ArgType5 , class ArgType6 , class ArgType7 >
KOKKOS_INLINE_FUNCTION
DstViewType
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
  DstViewType dst ;

  Impl::ViewAssignment<typename DstViewType::specialize,S>( dst, src, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7 );

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
#if KOKKOS_MIRROR_VIEW_OPTIMIZE
  return create_mirror_view( view );
#else
  typedef typename View<T,L,D,M,S>::HostMirror host_view ;
  host_view tmp ;
  Impl::ViewAssignment< S >( tmp , view );
  return tmp ;
#endif
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/Kokkos_ViewDefault.hpp>
#include <impl/Kokkos_ViewScalar.hpp>
#include <impl/Kokkos_ViewTileLeft.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

