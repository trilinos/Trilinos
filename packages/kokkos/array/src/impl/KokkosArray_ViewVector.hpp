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

#ifndef KOKKOSARRAY_VIEWVECTOR_HPP
#define KOKKOSARRAY_VIEWVECTOR_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

struct LayoutVector {};

template< class ViewTraits , class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits , ValueType , LayoutLeft , 1 , MemorySpace , MemoryTraits , void >
{ typedef LayoutVector type ; };

template< class ViewTraits , class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ViewTraits , ValueType , LayoutRight , 1 , MemorySpace , MemoryTraits , void >
{ typedef LayoutVector type ; };

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutVector , void , void >
{
  template< class T , class L , class D , class M >
  KOKKOSARRAY_INLINE_FUNCTION static
  size_t allocation_count( const View<T,L,D,M,LayoutVector> & dst )
  { return dst.m_shape.N0 ; }


  template< class T , class L , class D , class M >
  ViewAssignment( View<T,L,D,M,LayoutVector> & dst ,
                  const typename enable_if< ViewTraits<T,L,D,M>::is_managed , std::string >::type & label ,
                  const size_t n0 )
  {
    typedef View<T,L,D,M,LayoutVector> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::shape_type    shape_type ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape , n0 );

    dst.m_ptr_on_device = (typename DstViewType::value_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::value_type) ,
                              sizeof(typename DstViewType::value_type) ,
                              dst.m_shape.N0 );

    ViewInitialize< DstViewType >::apply( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  ViewAssignment(       View<DT,DL,DD,DM,LayoutVector> & dst ,
                  const View<ST,SL,SD,SM,LayoutVector> & src ,
                  typename enable_if<
                    is_same< View<DT,DL,DD,DM,LayoutVector> ,
                             typename View<ST,SL,SD,SM,LayoutVector>::HostMirror >::value
                  >::type * = 0 )
  {
    ViewAssignment( dst , "mirror" , src.m_shape.N0 );
  }
};

template<>
struct ViewAssignment< LayoutVector , LayoutVector , void >
{
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutVector> & dst , 
                  const View<ST,SL,SD,SM,LayoutVector> & src ,
                  typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                  ) >::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;
    typedef typename traits_type::shape_type shape_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape , src.m_shape.N0 );
    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }


  /** \brief  Vector subrange of a vector */
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM ,
            typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutVector> & dst , 
                  const View<ST,SL,SD,SM,LayoutVector> & src ,
                  const std::pair<iType,iType> & range ,
                  typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                    &&
                    ( ViewTraits<DT,DL,DD,DM>::rank_dynamic == 1 )
                  ) >::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;
    typedef typename traits_type::shape_type shape_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0      = 0 ;
    dst.m_ptr_on_device = 0 ;

    if ( range.first < range.second ) {
      assert_shape_bounds( src.m_shape , range.first );
      assert_shape_bounds( src.m_shape , range.second - 1 );

      dst.m_shape.N0 = range.second - range.first ;
      dst.m_ptr_on_device = src.m_ptr_on_device + range.first ;

      ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
    }
  }
};

template<>
struct ViewAssignment< LayoutScalar , LayoutVector , void >
{
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst , 
                  const View<ST,SL,SD,SM,LayoutVector> & src ,
                  const typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                  ), unsigned >::type i0 )
  {
    typedef View<DT,DL,DD,DM> traits_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }
};

template<>
struct ViewAssignment< LayoutVector , LayoutLeft , void >
{
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutVector> & dst , 
                  const View<ST,SL,SD,SM,LayoutLeft> & src ,
                  const typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                  ), unsigned >::type i1 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_shape.N0      = src.m_shape.N0 ;
    dst.m_ptr_on_device = src.m_ptr_on_device + src.m_stride * i1 ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType , class DeviceType , class MemoryTraits >
class View< DataType , LayoutType , DeviceType , MemoryTraits , Impl::LayoutVector >
  : public ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits >
{
private:

  template< class , class , class > friend class Impl::ViewAssignment ;

  typedef ViewTraits< DataType , LayoutType , DeviceType , MemoryTraits > traits ;

  typedef Impl::ViewAssignment<Impl::LayoutVector> alloc ;
  typedef Impl::ViewAssignment<Impl::LayoutVector,Impl::LayoutVector> assign ;

  typename traits::value_type * m_ptr_on_device ;
  typename traits::shape_type   m_shape ;

public:

  typedef Impl::LayoutVector specialize ;

  typedef View< typename traits::const_data_type ,
                typename traits::layout_type ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::layout_type ,
                Host > HostMirror ;

  enum { Rank = 1 };

  KOKKOSARRAY_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_1() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_2() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_3() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_4() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_5() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_6() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_7() const { return 1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  KOKKOSARRAY_INLINE_FUNCTION
  View() : m_ptr_on_device(0) {}

  KOKKOSARRAY_INLINE_FUNCTION
  ~View() { Impl::ViewTracking< traits >::decrement( m_ptr_on_device ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View( const View & rhs )
    : m_ptr_on_device(0) { assign( *this , rhs ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { assign( *this , rhs ); return *this ; }

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,Impl::LayoutVector> & rhs )
    : m_ptr_on_device(0) { assign( *this , rhs ); }

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,Impl::LayoutVector> & rhs )
    { assign( *this , rhs ); return *this ; }

  //------------------------------------

  explicit
  View( const std::string & label , const unsigned n0 = 0 ) : m_ptr_on_device(0)
    { alloc( *this , label , n0 ); }

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<( 0 == iType(0) ), typename traits::value_type >::type &
  operator()( const iType & i0 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[i0] ;
    }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  typename Impl::enable_if<( 0 == iType(0) ), typename traits::value_type >::type &
  operator[]( const iType & i0 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape, i0 );
      KOKKOSARRAY_ASSUME_ALIGNED( typename traits::memory_space , m_ptr_on_device );

      return m_ptr_on_device[i0] ;
    }
};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWVECTOR_HPP */

