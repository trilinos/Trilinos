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

#ifndef KOKKOSARRAY_ARRAYVIEWOPERRIGHT_HPP
#define KOKKOSARRAY_ARRAYVIEWOPERRIGHT_HPP

#include <KokkosArray_Macros.hpp>
#include <KokkosArray_View.hpp>
#include <impl/KokkosArray_Shape.hpp>

namespace KokkosArray {
namespace Impl {

struct LayoutArrayRight {};

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 1 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 2 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 3 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 4 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 5 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 6 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 7 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

template< class ViewTraits , typename ScalarType , unsigned Count ,
          class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits ,
                       Array< ScalarType , Count , void > ,
                       LayoutRight , 8 , MemorySpace , MemoryTraits ,
                       void >
{ typedef LayoutArrayRight type ; };

//-----------------------------------------------------------------------------

template<>
struct ViewAssignment<LayoutArrayRight,void,void>
{
  typedef LayoutArrayRight Specialize ;

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  size_t allocation_count( const View<T,L,D,M,Specialize> & dst )
  {
    return dst.m_shape.N0 * dst.m_stride ;
  }

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  void decrement( View<T,L,D,M,Specialize> & dst )
  {
    typedef View<T,L,D,M,Specialize> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< memory_space , memory_traits >::decrement( dst.m_ptr_on_device );
  }

  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  void increment( View<T,L,D,M,Specialize> & dst )
  {
    typedef View<T,L,D,M,Specialize> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< memory_space , memory_traits >::increment( dst.m_ptr_on_device );
  }

private:

  template< class T , class L , class D , class M >
  inline
  void allocate( View<T,L,D,M,Specialize> & dst , const std::string & label )
  {
    typedef View<T,L,D,M,Specialize> DstViewType ;
    typedef typename DstViewType::scalar_type   scalar_type ;
    typedef typename DstViewType::memory_space  memory_space ;

    const size_t count = allocation_count( dst );

    dst.m_ptr_on_device = (scalar_type *)
      memory_space::allocate( label , typeid(scalar_type) , sizeof(scalar_type) , count );

    ViewInitialize< DstViewType >::apply( dst );
  }

public:

  template< class T , class L , class D , class M >
  inline
  ViewAssignment( View<T,L,D,M,Specialize> & dst ,
                  const typename enable_if< ViewTraits<T,L,D,M>::is_managed , std::string >::type & label ,
                  const size_t n0 = 0 ,
                  const size_t n1 = 0 ,
                  const size_t n2 = 0 ,
                  const size_t n3 = 0 ,
                  const size_t n4 = 0 ,
                  const size_t n5 = 0 ,
                  const size_t n6 = 0 ,
                  const size_t n7 = 0 )
  {
    typedef View<T,L,D,M,Specialize> DstViewType ;
    typedef typename DstViewType::scalar_type   scalar_type ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;

    decrement( dst );

    shape_type::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );

    dst.m_stride =
      memory_space::preferred_alignment( dst.m_shape.scalar_size ,
                                         dst.m_shape.N1 * dst.m_shape.N2 *
                                         dst.m_shape.N3 * dst.m_shape.N4 *
                                         dst.m_shape.N5 * dst.m_shape.N6 *
                                         dst.m_shape.N7 );

    allocate( dst , label );
  }

  template< class T , class L , class D , class M >
  inline
  ViewAssignment( View<T,L,D,M,Specialize> & dst ,
                  const typename enable_if< ViewTraits<T,L,D,M>::is_managed , std::string >::type & label ,
                  const typename ViewTraits<T,L,D,M>::shape_type shape )
  {
    ViewAssignment( dst, label, shape.N0, shape.N1, shape.N2, shape.N3,
                                shape.N4, shape.N5, shape.N6 );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  typename enable_if<
                    is_same< View<DT,DL,DD,DM,Specialize> ,
                             typename View<ST,SL,SD,SM,Specialize>::HostMirror >::value
                  >::type * = 0 )
  {
    dst.m_shape  = src.m_shape ;
    dst.m_stride = src.m_stride ;
    allocate( dst , "mirror" );
  }
};

template<>
struct ViewAssignment< LayoutArrayRight , LayoutArrayRight , void >
{
  typedef LayoutArrayRight Specialize ;

  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Specialize> & dst ,
                  const View<ST,SL,SD,SM,Specialize> & src ,
                  const typename enable_if<(
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ShapeCompatible< typename ViewTraits<DT,DL,DD,DM>::shape_type ,
                                     typename ViewTraits<ST,SL,SD,SM>::shape_type >::value
                  )>::type * = 0 )
  {
    typedef View<DT,DL,DD,DM,Specialize> DstViewType ;
    typedef typename DstViewType::shape_type    shape_type ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewAssignment< Specialize >::decrement( dst );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 );

    dst.m_stride        = src.m_stride ;
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewAssignment< Specialize >::increment( dst );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType , class DeviceType , class MemoryTraits >
class View< DataType , LayoutType , DeviceType , MemoryTraits , Impl::LayoutArrayRight >
  : public ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits >
{
private:

  template< class , class , class > friend class Impl::ViewAssignment ;

  typedef ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits > traits ;

  typedef Impl::ViewAssignment<Impl::LayoutArrayRight> alloc ;
  typedef Impl::ViewAssignment<Impl::LayoutArrayRight,Impl::LayoutArrayRight> assign ;

  typedef Array< typename traits::scalar_type ,
                 ArrayTypeCount< typename traits::value_type >::value ,
                 ArrayProxyContiguous > ArrayProxyType ;

  typename traits::scalar_type * m_ptr_on_device ;
  unsigned                       m_stride ;
  typename traits::shape_type    m_shape ;

public:

  typedef Impl::LayoutArrayRight specialize ;

  typedef View< typename traits::const_data_type ,
                typename traits::layout_type ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::layout_type ,
                Host > HostMirror ;

  enum { Rank = traits::rank };

  KOKKOSARRAY_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_shape.N2 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_shape.N3 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_shape.N4 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_shape.N5 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_shape.N6 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_shape.N7 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  KOKKOSARRAY_INLINE_FUNCTION
  View() : m_ptr_on_device(0), m_stride(0) { traits::shape_type::assign(m_shape,0,0,0,0,0,0,0,0); }

  KOKKOSARRAY_INLINE_FUNCTION
  ~View() { alloc::decrement( *this ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View( const View & rhs ) : m_ptr_on_device(0) { assign( *this , rhs ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { assign( *this , rhs ); return *this ; }

  //------------------------------------

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,specialize> & rhs )
    : m_ptr_on_device(0) { assign( *this , rhs ); }

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,specialize> & rhs )
    { assign( *this , rhs ); return *this ; }

  //------------------------------------

  explicit
  View( const std::string & label ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    { alloc( *this, label, n0, n1, n2, n3, n4, n5, n6, n7 ); }

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::scalar_type * ptr_on_device() const { return m_ptr_on_device ; }

  //------------------------------------

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) &&
      1 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i0 ];
    }
  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) &&
      2 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,i1 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i1 + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) &&
      3 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,i2 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i2 + m_shape.N2 * (
                              i1 ) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      4 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,i3 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_stride ];
    }
  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      0 == iType4(0) &&
      5 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,i4 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      0 == iType4(0) && 0 == iType5(0) &&
      6 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,i5 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      0 == iType4(0) && 0 == iType5(0) && 0 == iType6(0) &&
      7 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      0 == iType4(0) && 0 == iType5(0) && 0 == iType6(0) && 0 == iType7(0) &&
      8 == traits::rank
    ), typename traits::scalar_type >::type & operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[ i7 + m_shape.N7 * (
                              i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))))) + i0 * m_stride ];
    }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator * () const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device , 1 );
    }

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) &&
      2 == traits::rank
    ), ArrayProxyType >::type operator ()
    ( const iType0 & i0 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape, i0,0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device + m_stride * i0 );
    }

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) &&
      3 == traits::rank
    ), ArrayProxyType >::type operator ()
    ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( m_shape, i0,i1,0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N2 * ( i1 ) +
                             m_stride * i0 );
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) &&
      4 == traits::rank
    ), ArrayProxyType >::type operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( m_shape, i0,i1,i2,0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 )) +
                             m_stride * i0 );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      5 == traits::rank
    ), ArrayProxyType >::type operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( m_shape, i0,i1,i2,i3,0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 ))) +
                             m_stride * i0 );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      0 == iType4(0) &&
      6 == traits::rank
    ), ArrayProxyType >::type operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( m_shape, i0,i1,i2,i3,i4,0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N5 * ( i4 +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 )))) +
                             m_stride * i0 );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      0 == iType4(0) && 0 == iType5(0) &&
      7 == traits::rank
    ), ArrayProxyType >::type operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( m_shape, i0,i1,i2,i3,i4,i5,0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N6 * ( i5 +
                             m_shape.N5 * ( i4 +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 ))))) +
                             m_stride * i0 );
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<(
      0 == iType0(0) && 0 == iType1(0) && 0 == iType2(0) && 0 == iType3(0) &&
      0 == iType4(0) && 0 == iType5(0) && 0 == iType6(0) &&
      8 == traits::rank
    ), ArrayProxyType >::type operator ()
    ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
      const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( m_shape, i0,i1,i2,i3,i4,i5,i6,0 );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N7 * ( i6 +
                             m_shape.N6 * ( i5 +
                             m_shape.N5 * ( i4 +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 )))))) +
                             m_stride * i0 );
    }
};

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_ARRAYVIEWOPERLEFT_HPP */

