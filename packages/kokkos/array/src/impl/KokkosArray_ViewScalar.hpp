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

#ifndef KOKKOSARRAY_VIEWSCALAR_HPP
#define KOKKOSARRAY_VIEWSCALAR_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

/** \brief  Specialization for a Rank = 0 array */
struct LayoutScalar ;

template< typename ScalarType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ScalarType , ScalarType ,
                       LayoutLeft , unsigned_<0> , unsigned_<0> ,
                       MemorySpace , MemoryTraits >
{ typedef LayoutScalar type ; };

template< typename ScalarType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ScalarType , ScalarType ,
                       LayoutRight , unsigned_<0> , unsigned_<0> ,
                       MemorySpace , MemoryTraits >
{ typedef LayoutScalar type ; };

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< LayoutScalar , void , void >
{
  template< class T , class L , class D , class M >
  ViewAssignment( View<T,L,D,M,LayoutScalar> & dst ,
                  typename enable_if< (
                    is_same< typename ViewTraits<T,L,D,M>::memory_traits ,
                             MemoryManaged >::value
                  ) , const std::string >::type & label )
  {
    typedef View<T,L,D,M,LayoutScalar> DstViewType ;
    typedef typename DstViewType::memory_space  memory_space ;
    typedef typename DstViewType::memory_traits memory_traits ;

    ViewTracking< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = (typename DstViewType::value_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::value_type) ,
                              sizeof(typename DstViewType::value_type) ,
                              1 );

    ViewInitialize< DstViewType >::apply( dst );
  }

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutScalar> & ,
                  typename enable_if< (
                    is_same< View<DT,DL,DD,DM,LayoutScalar> ,
                             typename View<ST,SL,SD,SM,LayoutScalar>::HostMirror >::value
                  ) >::type * = 0 )
  {
    ViewAssignment( dst , "mirror" );
  }
};

template<>
struct ViewAssignment< LayoutScalar , LayoutScalar , void >
{
  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( View<DT,DL,DD,DM,LayoutScalar> & dst , 
                  const View<ST,SL,SD,SM,LayoutScalar> & src ,
                  typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value
                  ) >::type * = 0 )
  {
    typedef ViewTraits<DT,DL,DD,DM> traits ;
    typedef typename traits::memory_space  memory_space ;
    typedef typename traits::memory_traits memory_traits ;

    ViewTracking< traits >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    ViewTracking< traits >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------
// Assign scalar from default layout array

template<>
struct ViewAssignment< LayoutScalar , LayoutDefault , void >
{
  //------------------------------------
  // Scalar from Rank=1

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const typename enable_if< (
                    ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                     ViewTraits<ST,SL,SD,SM> >::value &&
                    ( ViewTraits<ST,SL,SD,SM>::rank == 1 )
                  ), unsigned >::type i0 )
  {
    typedef View<DT,DL,DD,DM> traits_type ;

    ViewTracking< traits_type >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 ;

    ViewTracking< traits_type >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=2

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 2 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i1 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape , i0 , i1 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device = src.m_ptr_on_device + i0 + src.m_stride * i1 ;
    }
    else {
      dst.m_ptr_on_device = src.m_ptr_on_device + i1 + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=3

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 3 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i2 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * i2 );
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
        i2 + src.m_shape.N2 * (
        i1 ) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=4

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar>  & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault> & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 4 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i3 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * i3 ));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
        i3 + src.m_shape.N3 * (
        i2 + src.m_shape.N2 * (
        i1 )) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=5

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 5 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i4 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * i4 )));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 ))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=6

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 6 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i5 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * i5 ))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 )))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=7

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 7 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i6 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5, i6 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * (
          i5 + src.m_shape.N5 * i6 )))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i6 + src.m_shape.N6 * (
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 ))))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }

  //------------------------------------
  // Scalar from Rank=8

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,LayoutScalar> & dst ,
                  const View<ST,SL,SD,SM,LayoutDefault>   & src ,
                  const unsigned i0 ,
                  const unsigned i1 ,
                  const unsigned i2 ,
                  const unsigned i3 ,
                  const unsigned i4 ,
                  const unsigned i5 ,
                  const unsigned i6 ,
                  const typename enable_if< (
                    ( ViewTraits<ST,SL,SD,SM>::rank == 8 )
                    &&
                    ( ValueCompatible< ViewTraits<DT,DL,DD,DM> ,
                                       ViewTraits<ST,SL,SD,SM> >::value )
                  ) , unsigned >::type i7 )
  {
    typedef ViewTraits<DT,DL,DD,DM> view_traits ;

    enum { is_left = is_same< typename ViewTraits<ST,SL,SD,SM>::array_layout , LayoutLeft >::value };

    assert_shape_bounds( src.m_shape, i0, i1, i2, i3, i4, i5, i6, i7 );

    ViewTracking< view_traits >::decrement( dst.m_ptr_on_device );

    if ( is_left ) {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i0 + src.m_stride * (
          i1 + src.m_shape.N1 * (
          i2 + src.m_shape.N2 * (
          i3 + src.m_shape.N3 * (
          i4 + src.m_shape.N4 * (
          i5 + src.m_shape.N5 * (
          i6 + src.m_shape.N6 * i7 ))))));
    }
    else {
      dst.m_ptr_on_device =
        src.m_ptr_on_device +
          i7 + src.m_shape.N7 * (
          i6 + src.m_shape.N6 * (
          i5 + src.m_shape.N5 * (
          i4 + src.m_shape.N4 * (
          i3 + src.m_shape.N3 * (
          i2 + src.m_shape.N2 * (
          i1 )))))) + i0 * src.m_stride ;
    }

    ViewTracking< view_traits >::increment( dst.m_ptr_on_device );
  }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class Arg1Type , class Arg2Type , class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Impl::LayoutScalar >
  : public ViewTraits< DataType , Arg1Type , Arg2Type , Arg3Type >
{
private:

  template< class , class , class > friend class Impl::ViewAssignment ;

  typedef ViewTraits< DataType , Arg1Type , Arg2Type , Arg3Type > traits ;

  typedef Impl::ViewAssignment< Impl::LayoutScalar > alloc ;
  typedef Impl::ViewAssignment< Impl::LayoutScalar ,
                                Impl::LayoutScalar > assign ;

  typename traits::value_type * m_ptr_on_device ;

public:

  typedef View< typename traits::const_data_type,
                typename traits::array_layout,
                typename traits::device_type,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                Host ,
                void > HostMirror ;

  enum { Rank = 0 };

  KOKKOSARRAY_INLINE_FUNCTION typename traits::shape_type shape() const { return typename traits::shape_type(); }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_0() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_1() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_2() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_3() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_4() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_5() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_6() const { return 1 ; }
  KOKKOSARRAY_INLINE_FUNCTION typename traits::size_type dimension_7() const { return 1 ; }

  template< typename iType >
  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & ) const { return 1 ; }

  KOKKOSARRAY_INLINE_FUNCTION
  View() : m_ptr_on_device(0) {}

  KOKKOSARRAY_INLINE_FUNCTION
  ~View()
    { Impl::ViewTracking< traits >::decrement( m_ptr_on_device ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View( const View & rhs )
    : m_ptr_on_device(0) { assign( *this , rhs ); }

  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View & rhs ) { assign( *this , rhs ); return *this ; }

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,Impl::LayoutScalar> & rhs )
    : m_ptr_on_device(0) { assign( *this , rhs ); }

  template< class RT , class RL , class RD , class RM >
  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,Impl::LayoutScalar> & rhs )
    { assign( *this , rhs ); return *this ; }

  //------------------------------------
  /** \brief  Allocation constructor */
  explicit
  View( const std::string & label ) : m_ptr_on_device(0) { alloc( *this , label ); }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  const View & operator = ( const typename traits::value_type & rhs ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      *m_ptr_on_device = rhs ;
      return *this ;
    }

  KOKKOSARRAY_INLINE_FUNCTION
  View & operator = ( const typename traits::value_type & rhs )
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      *m_ptr_on_device = rhs ;
      return *this ;
    }

  KOKKOSARRAY_INLINE_FUNCTION
  operator typename traits::value_type & () const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      return *m_ptr_on_device ;
    }

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::value_type & operator()() const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      return *m_ptr_on_device ;
    }

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::value_type & operator*() const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );
      return *m_ptr_on_device ;
    }

  //------------------------------------

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  KOKKOSARRAY_INLINE_FUNCTION
  typename traits::size_type capacity() const { return 1 ; }
};

} /* namespace KokkosArray */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_VIEWSCALAR_HPP */

